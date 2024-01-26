#include "stencil/boundary_matrix.hpp"
#include "stencil_slave.hpp"

#include <algorithm>

namespace {
template <class ElemTy>
class SlaveBlock : public BoundaryMatrixView<ElemTy> {
    using Base = BoundaryMatrixView<ElemTy>;

    enum BoundaryPosition {
        Left,
        Right,
#ifdef PACK_ROW
        Top,
        Bottom,
#endif
        BoundarySize,
    };

public:
    /// Construct the slave block from the host block.
    explicit SlaveBlock(BoundaryMatrixView<ElemTy> const& host_block) :
        Base(
            /*width=*/host_block.width(),
            /*height=*/host_block.height(),
            /*neighbor_width=*/host_block.boundary_width(),
            /*neighbor_height=*/host_block.boundary_height(),
            // We must set the stride of the new block as the actual width of the origin block,
            // because we will allocate contiguous memory for the new block.
            /*stride=*/host_block.width_with_boundary(),
            // We cannot set data here because we want to use the `elem_with_boundary_count()`
            // method in the base class to compute the number of elements, which requires the base
            // object to be initialized.
            /*data=*/nullptr
        ),
        temp_buf() {
        // Modify the data member in the base class.
        Base::_data = (ElemTy*) ldm_malloc(Base::elem_with_boundary_count() * sizeof(ElemTy));

        std::size_t const temp_buf_size = Base::boundary_width() * Base::height() * sizeof(ElemTy);
        std::generate_n(std::begin(temp_buf), 2, [&] {
            return (ElemTy*) ldm_malloc(temp_buf_size);
        });

#ifdef PACK_ROW
        {
            std::size_t const temp_buf_size =
                Base::boundary_height() * Base::width() * sizeof(ElemTy);
            std::generate_n(&temp_buf[Top], 2, [&] {
                return (ElemTy*) ldm_malloc(temp_buf_size);
            });
        }
#endif  // PACK_ROW
    }

    // Since the class itself is actually an owner, we should provide appropriate move constructor
    // and move assignment. But I'm too lazy to implement them.
    SlaveBlock(SlaveBlock&&) = delete;
    auto operator=(SlaveBlock&&) -> SlaveBlock& = delete;

    ~SlaveBlock() {
        if (Base::_data)
            ldm_free(Base::_data, Base::elem_with_boundary_count() * sizeof(ElemTy));

        std::size_t const temp_buf_size = Base::boundary_width() * Base::height() * sizeof(ElemTy);
        for (unsigned i = 0; i != 2; ++i) {
            if (temp_buf[i]) {
                ldm_free(temp_buf[i], temp_buf_size);
            }
        }

#ifdef PACK_ROW
        {
            std::size_t const temp_buf_size =
                Base::boundary_height() * Base::width() * sizeof(ElemTy);
            std::for_each_n(&temp_buf[Top], 2, [&](ElemTy* arr) {
                if (arr) {
                    ldm_free(arr, temp_buf_size);
                }
            });
        }
#endif  // PACK_ROW
    }

    auto empty() const -> bool {
        return Base::_data == nullptr
            || std::any_of(std::begin(temp_buf), std::end(temp_buf), [](ElemTy* ptr) {
                   return ptr == nullptr;
               });
    }

    /// Loads all elements (including boundaries) to LDM.
    void dma_get_all_from(BoundaryMatrixView<ElemTy>& host_block) {
        athread_dma_get_stride(
            &Base::elem_with_boundary_at(0, 0),
            &host_block.elem_with_boundary_at(0, 0),
            Base::elem_with_boundary_count() * sizeof(ElemTy),
            Base::width_with_boundary() * sizeof(ElemTy),
            (host_block.internal_data_stride() - Base::width_with_boundary()) * sizeof(ElemTy)
        );
    }

    void dma_iget_top_bottom_boundaries_from(
        BoundaryMatrixView<ElemTy>& host_block,
        athread_rply_t* reply
    ) {
        std::size_t const boundary_elem_count =
            Base::boundary_height() * Base::width_with_boundary();
        // top
        athread_dma_iget_stride(
            &Base::elem_with_boundary_at(0, 0),
            &host_block.elem_with_boundary_at(0, 0),
            boundary_elem_count * sizeof(ElemTy),
            Base::width_with_boundary() * sizeof(ElemTy),
            (host_block.internal_data_stride() - Base::width_with_boundary()) * sizeof(ElemTy),
            reply
        );

        // bottom
        athread_dma_iget_stride(
            &Base::elem_with_boundary_at(
                Base::height_with_boundary() - Base::boundary_height(),
                0
            ),
            &host_block.elem_with_boundary_at(  //
                Base::height_with_boundary() - Base::boundary_height(),
                0
            ),
            boundary_elem_count * sizeof(ElemTy),
            Base::width_with_boundary() * sizeof(ElemTy),
            (host_block.internal_data_stride() - Base::width_with_boundary()) * sizeof(ElemTy),
            reply
        );
    }

    /// Reads the left and right boundaries into the temp buffers by DMA.
    void dma_iget_left_right_boundaries_to_temp_buffer_from(
        BoundaryMatrixView<ElemTy>& host_block,
        athread_rply_t* reply
    ) {
        std::size_t const boundary_elem_count = Base::boundary_width() * Base::height();

        // Left
        athread_dma_iget_stride(
            temp_buf[0],
            &host_block.elem_with_boundary_at(Base::boundary_height(), 0),
            boundary_elem_count * sizeof(ElemTy),
            Base::boundary_width() * sizeof(ElemTy),
            (host_block.internal_data_stride() - Base::boundary_width()) * sizeof(ElemTy),
            reply
        );

        // Right
        athread_dma_iget_stride(
            temp_buf[1],
            &host_block.elem_with_boundary_at(
                Base::boundary_height(),
                Base::width_with_boundary() - Base::boundary_width()
            ),
            boundary_elem_count * sizeof(ElemTy),
            Base::boundary_width() * sizeof(ElemTy),
            (host_block.internal_data_stride() - Base::boundary_width()) * sizeof(ElemTy),
            reply
        );
    }

    /// Copies the content in the temp buffers to the left and right boundaries.
    void copy_temp_buf_to_boundaries() {
        auto* lbuf_ptr = temp_buf[0];
        auto* rbuf_ptr = temp_buf[1];
        std::size_t const row_elem_count = Base::boundary_width();
        std::size_t const right_boundary_col_beg =
            Base::width_with_boundary() - Base::boundary_width();

        for (unsigned row = Base::boundary_height(),
                      end = Base::height_with_boundary() - Base::boundary_height();
             row != end;
             ++row)
        {
            std::copy_n(lbuf_ptr, row_elem_count, &Base::elem_with_boundary_at(row, 0));
            std::copy_n(
                rbuf_ptr,
                row_elem_count,
                &Base::elem_with_boundary_at(row, right_boundary_col_beg)
            );

            lbuf_ptr += row_elem_count;
            rbuf_ptr += row_elem_count;
        }
    }

    auto dma_iput_top_boundary_to(  //
        BoundaryMatrixView<ElemTy>& host_block,
        athread_rply_t* reply
    ) -> unsigned {
#ifdef PACK_ROW
        // TODO: Implementation.
#else
        if (Base::boundary_height() == 1) {
            // Fast path
            athread_dma_iput(
                &host_block.elem_at(0, 0),
                &Base::elem_at(0, 0),
                Base::width() * sizeof(ElemTy),
                reply
            );
        } else {
            for (unsigned r = 0; r != Base::boundary_height(); ++r) {
                athread_dma_iput(
                    &host_block.elem_at(r, 0),
                    &Base::elem_at(r, 0),
                    Base::width() * sizeof(ElemTy),
                    reply
                );
            }
        }

        return Base::boundary_height();
#endif  // PACK_ROW
    }

    auto dma_iput_bottom_boundary_to(  //
        BoundaryMatrixView<ElemTy>& host_block,
        athread_rply_t* reply
    ) -> unsigned {
#ifdef PACK_ROW
        // TODO: Implementation.
#else
        if (Base::boundary_height() == 1) {
            // Fast path
            athread_dma_iput(
                &host_block.elem_at(Base::height() - 1, 0),
                &Base::elem_at(Base::height() - 1, 0),
                Base::width() * sizeof(ElemTy),
                reply
            );
        } else {
            for (unsigned r = 0; r != Base::boundary_height(); ++r) {
                athread_dma_iput(
                    &host_block.elem_at(Base::height() - 1 - r, 0),
                    &Base::elem_at(Base::height() - 1 - r, 0),
                    Base::width() * sizeof(ElemTy),
                    reply
                );
            }
        }

        return Base::boundary_height();
#endif  // PACK_ROW
    }

    void dma_iput_left_boundary_to(BoundaryMatrixView<ElemTy>& host_block, athread_rply_t* reply) {
        // Copy the boundaries to temp buffers.
        auto* lbuf_ptr = temp_buf[0];
        std::size_t const row_elem_count = Base::boundary_width();

        for (unsigned row = Base::boundary_height(),
                      end = Base::height() - Base::boundary_height();
             row != end;
             ++row)
        {
            std::copy_n(&Base::elem_at(row, 0), row_elem_count, lbuf_ptr);
            lbuf_ptr += row_elem_count;
        }

        // DMA to the host memory.
        athread_dma_iput_stride(
            &host_block.elem_at(Base::boundary_height(), 0),
            temp_buf[0],
            (lbuf_ptr - temp_buf[0]) * sizeof(ElemTy),
            Base::boundary_width() * sizeof(ElemTy),
            (host_block.internal_data_stride() - Base::boundary_width()) * sizeof(ElemTy),
            reply
        );
    }

    void dma_iput_right_boundary_to(
        BoundaryMatrixView<ElemTy>& host_block,
        athread_rply_t* reply
    ) {
        // Copy the boundaries to temp buffers.
        auto* rbuf_ptr = temp_buf[1];
        std::size_t const row_elem_count = Base::boundary_width();
        std::size_t const right_boundary_col_beg = Base::width() - Base::boundary_width();

        for (unsigned row = Base::boundary_height(),
                      end = Base::height() - Base::boundary_height();
             row != end;
             ++row)
        {
            std::copy_n(&Base::elem_at(row, right_boundary_col_beg), row_elem_count, rbuf_ptr);
            rbuf_ptr += row_elem_count;
        }

        // DMA to the host memory.
        athread_dma_iput_stride(
            &host_block.elem_at(Base::boundary_height(), right_boundary_col_beg),
            temp_buf[1],
            (rbuf_ptr - temp_buf[1]) * sizeof(ElemTy),
            Base::boundary_width() * sizeof(ElemTy),
            (host_block.internal_data_stride() - Base::boundary_width()) * sizeof(ElemTy),
            reply
        );
    }

    void dma_put_content_to(BoundaryMatrixView<ElemTy>& host_block) {
        for (unsigned r = 0; r != Base::height(); ++r) {
            athread_dma_put(
                &host_block.elem_at(r, 0),
                &Base::elem_at(r, 0),
                Base::width() * sizeof(ElemTy)
            );
        }
    }

    /// Swap two `SlaveBlock`s.
    ///
    /// Note that since there is no move constructor and move assignment in `SlaveBlock`, we
    /// provide this simple swap method. But this method assumes that the two objects being swapped
    /// are exactly the same size and shape.
    void swap(SlaveBlock& other) noexcept {
        std::swap(Base::_data, other._data);
    }

private:
    // Temp buffers used to store left and right boundaries to be sent or received.
    ElemTy* temp_buf[BoundarySize];
};
}  // namespace

void stencil_iterate_dma_slave_pack(Arguments* _args) {
    Arguments args;
    athread_dma_get(&args, _args, sizeof(Arguments));

    // Divide the host matrix into blocks.
    BoundaryMatrixView<float> local_host_input =
        args.input.block_subview(args.block_size, args.block_size, _ROW, _COL);
    BoundaryMatrixView<float> local_host_output =
        args.output.block_subview(args.block_size, args.block_size, _ROW, _COL);

    // Allocate LDM.
    SlaveBlock slave_input(local_host_input);
    SlaveBlock slave_output(local_host_output);

    if (slave_input.empty() || slave_output.empty()) {
        return;
    }

    // Load content into LDM.
    slave_input.dma_get_all_from(local_host_input);

    athread_rply_t dma_reply = 0;
    unsigned dma_reply_count = 0;

    // Iterations.
    for (unsigned i = 0; i != args.iterations; ++i) {
        if (i != 0) {
            // Load Boundaries.
            slave_input.dma_iget_left_right_boundaries_to_temp_buffer_from(  //
                local_host_input,
                &dma_reply
            );
            slave_input.dma_iget_top_bottom_boundaries_from(local_host_input, &dma_reply);

            dma_reply_count += 4;
        }

        // Compute contents.
        for (unsigned r = slave_input.boundary_height(),
                      r_end = slave_input.height() - slave_input.boundary_height();
             r != r_end;
             ++r)
        {
            for (unsigned c = slave_input.boundary_width(),
                          c_end = slave_input.width() - slave_input.boundary_width();
                 c != c_end;
                 ++c)
            {
                slave_output.elem_at(r, c) = 0.25f
                    * (slave_input.elem_at(r - 1, c) + slave_input.elem_at(r, c - 1)
                       + slave_input.elem_at(r, c + 1) + slave_input.elem_at(r + 1, c));
            }
        }

        if (i != 0) {
            // Wait for boundaries.
            athread_dma_wait_value(&dma_reply, dma_reply_count);
            // Copy the temp buffers to boundaries.
            slave_input.copy_temp_buf_to_boundaries();
        }

        // Top boundary.
        for (unsigned c = slave_input.boundary_width();
             c != slave_input.width_with_boundary() - slave_input.boundary_width();
             ++c)
        {
            unsigned const r = slave_input.boundary_height();
            slave_output.elem_with_boundary_at(r, c) = 0.25f
                * (slave_input.elem_with_boundary_at(r - 1, c)
                   + slave_input.elem_with_boundary_at(r, c - 1)
                   + slave_input.elem_with_boundary_at(r, c + 1)
                   + slave_input.elem_with_boundary_at(r + 1, c));
        }

        dma_reply_count += slave_output.dma_iput_top_boundary_to(local_host_output, &dma_reply);

        // Bottom boundary.
        for (unsigned c = slave_input.boundary_width();
             c != slave_input.width_with_boundary() - slave_input.boundary_width();
             ++c)
        {
            unsigned const r =
                slave_input.height_with_boundary() - slave_input.boundary_height() - 1;
            slave_output.elem_with_boundary_at(r, c) = 0.25f
                * (slave_input.elem_with_boundary_at(r - 1, c)
                   + slave_input.elem_with_boundary_at(r, c - 1)
                   + slave_input.elem_with_boundary_at(r, c + 1)
                   + slave_input.elem_with_boundary_at(r + 1, c));
        }

        dma_reply_count += slave_output.dma_iput_bottom_boundary_to(local_host_output, &dma_reply);

        // Left boundary.
        for (unsigned r = slave_input.boundary_height();
             r != slave_input.height_with_boundary() - slave_input.boundary_height();
             ++r)
        {
            unsigned const c = slave_input.boundary_width();
            slave_output.elem_with_boundary_at(r, c) = 0.25f
                * (slave_input.elem_with_boundary_at(r - 1, c)
                   + slave_input.elem_with_boundary_at(r, c - 1)
                   + slave_input.elem_with_boundary_at(r, c + 1)
                   + slave_input.elem_with_boundary_at(r + 1, c));
        }

        slave_output.dma_iput_left_boundary_to(local_host_output, &dma_reply);
        ++dma_reply_count;

        // Right boundary.
        for (unsigned r = slave_input.boundary_height();
             r != slave_input.height_with_boundary() - slave_input.boundary_height();
             ++r)
        {
            unsigned const c =
                slave_input.width_with_boundary() - slave_input.boundary_width() - 1;
            slave_output.elem_with_boundary_at(r, c) = 0.25f
                * (slave_input.elem_with_boundary_at(r - 1, c)
                   + slave_input.elem_with_boundary_at(r, c - 1)
                   + slave_input.elem_with_boundary_at(r, c + 1)
                   + slave_input.elem_with_boundary_at(r + 1, c));
        }

        slave_output.dma_iput_right_boundary_to(local_host_output, &dma_reply);
        ++dma_reply_count;

        // Swap the input and result.
        slave_input.swap(slave_output);
        std::swap(local_host_input, local_host_output);

        athread_dma_wait_value(&dma_reply, dma_reply_count);

        // Synchronize.
        if (i != args.iterations - 1)
            athread_ssync_array();
    }

    // Write the final result.
    slave_input.dma_put_content_to(local_host_input);
}
