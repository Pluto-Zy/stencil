#include "stencil/boundary_matrix.hpp"
#include "stencil_slave.hpp"
#include <algorithm>
#include <iterator>

namespace {
template <class ElemTy>
class SlaveBlock {
    enum BoundaryPosition {
        Top,
        Bottom,
        LeftIn,
        Left,
        RightIn,
        Right,
        BoundaryCount,
    };

public:
    /// Construct the slave block from the host block.
    explicit SlaveBlock(BoundaryMatrixView<ElemTy> const& host_block) :
        _rows(host_block.height()),
        _cols(host_block.width()),
        _boundary_height(host_block.boundary_height()),
        _boundary_width(host_block.boundary_width()),
        _data((ElemTy*) ldm_malloc(inner_data_mem_size())) {
        std::generate(std::begin(boundaries), std::begin(boundaries) + 2, [this] {
            return (ElemTy*) ldm_malloc(boundary_mem_size(BoundaryPosition::Top));
        });
        std::generate(std::begin(boundaries) + 2, std::end(boundaries), [this] {
            return (ElemTy*) ldm_malloc(boundary_mem_size(BoundaryPosition::Left));
        });
    }

    SlaveBlock(SlaveBlock&& other) noexcept :
        _rows(other._rows),
        _cols(other._cols),
        _boundary_height(other._boundary_height),
        _boundary_width(other._boundary_width),
        _data(other._data) {
        std::copy(
            std::begin(other.boundaries),
            std::end(other.boundaries),
            std::begin(boundaries)
        );
        std::fill(std::begin(other.boundaries), std::end(other.boundaries), nullptr);
        other._data = nullptr;
        other._rows = other._cols = other._boundary_width = other._boundary_height = 0;
    }

    auto operator=(SlaveBlock&& other) noexcept -> SlaveBlock& {
        if (this != &other) {
            free_resource();

            _data = other._data;
            std::copy(
                std::begin(other.boundaries),
                std::end(other.boundaries),
                std::begin(boundaries)
            );
            _rows = other._rows;
            _cols = other._cols;
            _boundary_height = other._boundary_height;
            _boundary_width = other._boundary_width;

            other._data = nullptr;
            std::fill(std::begin(other.boundaries), std::end(other.boundaries), nullptr);
            other._rows = other._cols = other._boundary_width = other._boundary_height = 0;
        }

        return *this;
    }

    ~SlaveBlock() {
        free_resource();
    }

    auto rows() const -> unsigned {
        return _rows;
    }

    auto cols() const -> unsigned {
        return _cols;
    }

    auto empty() const -> bool {
        return _data == nullptr
            || std::any_of(std::begin(boundaries), std::end(boundaries), [](ElemTy* ptr) {
                   return ptr == nullptr;
               });
    }

    auto data_at(unsigned row, unsigned col) const -> ElemTy& {
        return _data[row * (_cols - _boundary_width * 2) + col - _boundary_width];
    }

    auto top_boundary_at(unsigned row, unsigned col) const -> ElemTy& {
        return boundaries[Top][row * _cols + col];
    }

    auto bottom_boundary_at(unsigned row, unsigned col) const -> ElemTy& {
        return boundaries[Bottom][row * _cols + col];
    }

    auto left_boundary_at(unsigned row, unsigned col) const -> ElemTy& {
        return boundaries[Left][row * _boundary_width + col];
    }

    auto left_inner_boundary_at(unsigned row, unsigned col) const -> ElemTy& {
        return boundaries[LeftIn][row * _boundary_width + col];
    }

    auto right_boundary_at(unsigned row, unsigned col) const -> ElemTy& {
        return boundaries[Right][row * _boundary_width + col];
    }

    auto right_inner_boundary_at(unsigned row, unsigned col) const -> ElemTy& {
        return boundaries[RightIn][row * _boundary_width + col];
    }

    void dma_get_content_from(BoundaryMatrixView<ElemTy>& src_block) const {
        unsigned const vec_block_size = _cols - 2 * _boundary_width;

        // main body
        athread_dma_get_stride(
            &data_at(0, _boundary_width),
            &src_block.elem_at(0, _boundary_width),
            inner_data_mem_size(),
            vec_block_size * sizeof(ElemTy),
            (src_block.internal_data_stride() - vec_block_size) * sizeof(ElemTy)
        );

        // inner left
        athread_dma_get_stride(
            &left_inner_boundary_at(0, 0),
            &src_block.elem_at(0, 0),
            boundary_mem_size(BoundaryPosition::LeftIn),
            _boundary_width * sizeof(ElemTy),
            (src_block.internal_data_stride() - _boundary_width) * sizeof(ElemTy)
        );

        // inner right
        athread_dma_get_stride(
            &right_inner_boundary_at(0, 0),
            &src_block.elem_at(0, src_block.width() - _boundary_width),
            boundary_mem_size(BoundaryPosition::RightIn),
            _boundary_width * sizeof(ElemTy),
            (src_block.internal_data_stride() - _boundary_width) * sizeof(ElemTy)
        );
    }

    void dma_iget_top_boundary_from(  //
        BoundaryMatrixView<ElemTy>& src_block,
        athread_rply_t* reply
    ) const {
        if (_boundary_height == 1) {
            athread_dma_iget(
                &top_boundary_at(0, 0),
                &src_block.elem_with_boundary_at(0, _boundary_width),
                boundary_mem_size(BoundaryPosition::Top),
                reply
            );
        } else {
            athread_dma_iget_stride(
                &top_boundary_at(0, 0),
                &src_block.elem_with_boundary_at(0, _boundary_width),
                boundary_mem_size(BoundaryPosition::Top),
                _cols * sizeof(ElemTy),
                (src_block.internal_data_stride() - _cols) * sizeof(ElemTy),
                reply
            );
        }
    }

    void dma_iget_bottom_boundary_from(  //
        BoundaryMatrixView<ElemTy>& src_block,
        athread_rply_t* reply
    ) const {
        if (_boundary_height == 1) {
            athread_dma_iget(
                &bottom_boundary_at(0, 0),
                &src_block.elem_with_boundary_at(
                    src_block.height_with_boundary() - src_block.boundary_height(),
                    _boundary_width
                ),
                boundary_mem_size(BoundaryPosition::Bottom),
                reply
            );
        } else {
            athread_dma_iget_stride(
                &bottom_boundary_at(0, 0),
                &src_block.elem_with_boundary_at(
                    src_block.height_with_boundary() - src_block.boundary_height(),
                    _boundary_width
                ),
                boundary_mem_size(BoundaryPosition::Bottom),
                _cols * sizeof(ElemTy),
                (src_block.internal_data_stride() - _cols) * sizeof(ElemTy),
                reply
            );
        }
    }

    void dma_iget_left_boundary_from(  //
        BoundaryMatrixView<ElemTy>& src_block,
        athread_rply_t* reply
    ) const {
        athread_dma_iget_stride(
            &left_boundary_at(0, 0),
            &src_block.elem_with_boundary_at(src_block.boundary_height(), 0),
            boundary_mem_size(BoundaryPosition::Left),
            _boundary_width * sizeof(ElemTy),
            (src_block.internal_data_stride() - _boundary_width) * sizeof(ElemTy),
            reply
        );
    }

    void dma_iget_right_boundary_from(  //
        BoundaryMatrixView<ElemTy>& src_block,
        athread_rply_t* reply
    ) const {
        athread_dma_iget_stride(
            &right_boundary_at(0, 0),
            &src_block.elem_with_boundary_at(
                src_block.boundary_height(),
                src_block.width_with_boundary() - src_block.boundary_width()
            ),
            boundary_mem_size(BoundaryPosition::Right),
            _boundary_width * sizeof(ElemTy),
            (src_block.internal_data_stride() - _boundary_width) * sizeof(ElemTy),
            reply
        );
    }

    void dma_iget_boundary_from(  //
        BoundaryMatrixView<ElemTy>& src_block,
        athread_rply_t* reply,
        unsigned reply_count
    ) const {
        dma_iget_top_boundary_from(src_block, reply);
        dma_iget_bottom_boundary_from(src_block, reply);
        dma_iget_left_boundary_from(src_block, reply);
        dma_iget_right_boundary_from(src_block, reply);

        reply_count += 4;
    }

    void dma_put_content_to(BoundaryMatrixView<ElemTy>& dst_block) const {
        unsigned const vec_block_size = _cols - 2 * _boundary_width;

        // main body
        athread_dma_put_stride(
            &dst_block.elem_at(0, _boundary_width),
            &data_at(0, _boundary_width),
            inner_data_mem_size(),
            vec_block_size * sizeof(ElemTy),
            (dst_block.internal_data_stride() - vec_block_size) * sizeof(ElemTy)
        );

        // inner left
        athread_dma_put_stride(
            &dst_block.elem_at(0, 0),
            &left_inner_boundary_at(0, 0),
            boundary_mem_size(BoundaryPosition::LeftIn),
            _boundary_width * sizeof(ElemTy),
            (dst_block.internal_data_stride() - _boundary_width) * sizeof(ElemTy)
        );

        // inner right
        athread_dma_put_stride(
            &dst_block.elem_at(0, dst_block.width() - _boundary_width),
            &right_inner_boundary_at(0, 0),
            boundary_mem_size(BoundaryPosition::RightIn),
            _boundary_width * sizeof(ElemTy),
            (dst_block.internal_data_stride() - _boundary_width) * sizeof(ElemTy)
        );
    }

    void dma_iput_top_boundary_to(  //
        BoundaryMatrixView<ElemTy>& dst_block,
        athread_rply_t* reply
    ) const {
        if (_boundary_height == 1) {
            athread_dma_iput(
                &dst_block.elem_at(0, _boundary_width),
                &data_at(0, _boundary_width),
                (_cols - _boundary_width * 2) * sizeof(ElemTy),
                reply
            );
        } else {
            unsigned const vec_block_size = _cols - _boundary_width * 2;
            athread_dma_iput_stride(
                &dst_block.elem_at(0, _boundary_width),
                &data_at(0, _boundary_width),
                vec_block_size * _boundary_height * sizeof(ElemTy),
                vec_block_size * sizeof(ElemTy),
                (dst_block.internal_data_stride() - vec_block_size) * sizeof(ElemTy),
                reply
            );
        }
    }

    void dma_iput_bottom_boundary_to(  //
        BoundaryMatrixView<ElemTy>& dst_block,
        athread_rply_t* reply
    ) const {
        if (_boundary_height == 1) {
            athread_dma_iput(
                &dst_block.elem_at(_rows - _boundary_height, _boundary_width),
                &data_at(_rows - _boundary_height, _boundary_width),
                (_cols - _boundary_width * 2) * sizeof(ElemTy),
                reply
            );
        } else {
            unsigned const vec_block_size = _cols - _boundary_width * 2;
            athread_dma_iput_stride(
                &dst_block.elem_at(_rows - _boundary_height, _boundary_width),
                &data_at(_rows - _boundary_height, _boundary_width),
                vec_block_size * _boundary_height * sizeof(ElemTy),
                vec_block_size * sizeof(ElemTy),
                (dst_block.internal_data_stride() - vec_block_size) * sizeof(ElemTy),
                reply
            );
        }
    }

    void dma_iput_left_boundary_to(  //
        BoundaryMatrixView<ElemTy>& dst_block,
        athread_rply_t* reply
    ) const {
        athread_dma_iput_stride(
            &dst_block.elem_at(0, 0),
            &left_inner_boundary_at(0, 0),
            boundary_mem_size(BoundaryPosition::LeftIn),
            _boundary_width * sizeof(ElemTy),
            (dst_block.internal_data_stride() - _boundary_width) * sizeof(ElemTy),
            reply
        );
    }

    void dma_iput_right_boundary_to(  //
        BoundaryMatrixView<ElemTy>& dst_block,
        athread_rply_t* reply
    ) const {
        athread_dma_iput_stride(
            &dst_block.elem_at(0, _cols - _boundary_width),
            &right_inner_boundary_at(0, 0),
            boundary_mem_size(BoundaryPosition::RightIn),
            _boundary_width * sizeof(ElemTy),
            (dst_block.internal_data_stride() - _boundary_width) * sizeof(ElemTy),
            reply
        );
    }

private:
    unsigned _rows;
    unsigned _cols;
    unsigned _boundary_height;
    unsigned _boundary_width;
    ElemTy* _data;
    ElemTy* boundaries[BoundaryCount];

    auto inner_data_mem_size() const -> std::size_t {
        return sizeof(ElemTy) * _rows * (_cols - _boundary_width * 2);
    }

    auto boundary_mem_size(BoundaryPosition position) const -> std::size_t {
        switch (position) {
        case Top:
        case Bottom:
            return _cols * _boundary_height * sizeof(ElemTy);
        case Left:
        case LeftIn:
        case Right:
        case RightIn:
            return _rows * _boundary_width * sizeof(ElemTy);
        default:
            __builtin_unreachable();
        }
    }

    void free_resource() noexcept {
        if (_data) {
            ldm_free(_data, inner_data_mem_size());
        }

        std::for_each(std::begin(boundaries), std::begin(boundaries) + 2, [this](ElemTy* ptr) {
            if (ptr) {
                ldm_free(ptr, boundary_mem_size(BoundaryPosition::Top));
            }
        });

        std::for_each(std::begin(boundaries) + 2, std::end(boundaries), [this](ElemTy* ptr) {
            if (ptr) {
                ldm_free(ptr, boundary_mem_size(BoundaryPosition::Left));
            }
        });
    }
};
}  // namespace

void stencil_iterate_dma(Arguments* _args) {
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
    slave_input.dma_get_content_from(local_host_input);

    athread_rply_t dma_reply = 0;
    unsigned dma_reply_count = 0;

    // Iterations.
    for (unsigned i = 0; i != args.iterations; ++i) {
        // Load boundaries.
        slave_input.dma_iget_boundary_from(local_host_input, &dma_reply, dma_reply_count);

        // Compute contents.
        for (unsigned row = 1; row != slave_input.rows() - 1; ++row) {
            unsigned col = 1;
            slave_output.data_at(row, col) = 0.25f
                * (slave_input.data_at(row - 1, col) + slave_input.left_inner_boundary_at(row, 0)
                   + slave_input.data_at(row, col + 1) + slave_input.data_at(row + 1, col));

            for (++col; col != slave_input.cols() - 2; ++col) {
                slave_output.data_at(row, col) = 0.25f
                    * (slave_input.data_at(row - 1, col) + slave_input.data_at(row, col - 1)
                       + slave_input.data_at(row, col + 1) + slave_input.data_at(row + 1, col));
            }

            slave_output.data_at(row, col) = 0.25f
                * (slave_input.data_at(row - 1, col) + slave_input.data_at(row, col - 1)
                   + slave_input.right_inner_boundary_at(row, 0)
                   + slave_input.data_at(row + 1, col));
        }

        // Wait for boundaries.
        athread_dma_wait_value(&dma_reply, dma_reply_count);

        {
            constexpr unsigned col = 0;
            unsigned row = 0;
            // Left-top corner.
            slave_output.left_inner_boundary_at(row, 0) = 0.25f
                * (slave_input.top_boundary_at(0, col) + slave_input.left_boundary_at(row, 0)
                   + slave_input.data_at(row, col + 1)
                   + slave_input.left_inner_boundary_at(row + 1, 0));

            // Left boundary.
            for (++row; row != slave_input.rows() - 1; ++row) {
                slave_output.left_inner_boundary_at(row, 0) = 0.25f
                    * (slave_input.left_inner_boundary_at(row - 1, 0)
                       + slave_input.left_boundary_at(row, 0) + slave_input.data_at(row, col + 1)
                       + slave_input.left_inner_boundary_at(row + 1, 0));
            }

            // Left-bottom corner.
            slave_output.left_inner_boundary_at(row, 0) = 0.25f
                * (slave_input.left_inner_boundary_at(row - 1, 0)
                   + slave_input.left_boundary_at(row, 0) + slave_input.data_at(row, col + 1)
                   + slave_input.bottom_boundary_at(0, col));
        }

        // Write left boundary to the main memory.
        slave_output.dma_iput_left_boundary_to(local_host_output, &dma_reply);

        {
            unsigned const col = slave_input.cols() - 1;
            unsigned row = 0;

            // Right-top corner.
            slave_output.right_inner_boundary_at(row, 0) = 0.25f
                * (slave_input.top_boundary_at(0, col) + slave_input.data_at(row, col - 1)
                   + slave_input.right_boundary_at(row, 0)
                   + slave_input.right_inner_boundary_at(row + 1, 0));

            // Right boundary.
            for (++row; row != slave_input.rows() - 1; ++row) {
                slave_output.right_inner_boundary_at(row, 0) = 0.25f
                    * (slave_input.right_inner_boundary_at(row - 1, 0)
                       + slave_input.data_at(row, col - 1) + slave_input.right_boundary_at(row, 0)
                       + slave_input.right_inner_boundary_at(row + 1, 0));
            }

            // Right-bottom corner.
            slave_output.right_inner_boundary_at(row, 0) = 0.25f
                * (slave_input.right_inner_boundary_at(row - 1, 0)
                   + slave_input.data_at(row, col - 1) + slave_input.right_boundary_at(row, 0)
                   + slave_input.bottom_boundary_at(0, col));
        }

        // Write right boundary to the main memory.
        slave_output.dma_iput_right_boundary_to(local_host_output, &dma_reply);

        // Top boundary.
        {
            unsigned col = 1;
            constexpr unsigned row = 0;

            slave_output.data_at(row, col) = 0.25f
                * (slave_input.top_boundary_at(0, col) + slave_input.left_inner_boundary_at(row, 0)
                   + slave_input.data_at(row, col + 1) + slave_input.data_at(row + 1, col));

            for (++col; col != slave_input.cols() - 2; ++col) {
                slave_output.data_at(row, col) = 0.25f
                    * (slave_input.top_boundary_at(0, col) + slave_input.data_at(row, col - 1)
                       + slave_input.data_at(row, col + 1) + slave_input.data_at(row + 1, col));
            }

            slave_output.data_at(row, col) = 0.25f
                * (slave_input.top_boundary_at(0, col) + slave_input.data_at(row, col - 1)
                   + slave_input.right_inner_boundary_at(row, 0)
                   + slave_input.data_at(row + 1, col));
        }

        // Write top boundary to the host memory.
        slave_output.dma_iput_top_boundary_to(local_host_output, &dma_reply);

        // Bottom boundary.
        {
            unsigned const row = slave_input.rows() - 1;
            unsigned col = 1;

            slave_output.data_at(row, col) = 0.25f
                * (slave_input.data_at(row - 1, col) + slave_input.left_inner_boundary_at(row, 0)
                   + slave_input.data_at(row, col + 1) + slave_input.bottom_boundary_at(0, col));

            for (++col; col != slave_input.cols() - 2; ++col) {
                slave_output.data_at(row, col) = 0.25f
                    * (slave_input.data_at(row - 1, col) + slave_input.data_at(row, col - 1)
                       + slave_input.data_at(row, col + 1)
                       + slave_input.bottom_boundary_at(0, col));
            }

            slave_output.data_at(row, col) = 0.25f
                * (slave_input.data_at(row - 1, col) + slave_input.data_at(row, col - 1)
                   + slave_input.right_inner_boundary_at(row, 0)
                   + slave_input.bottom_boundary_at(0, col));
        }

        // Write bottom boundary to the host memory.
        slave_output.dma_iput_bottom_boundary_to(local_host_output, &dma_reply);
        dma_reply_count += 4;

        // Swap the input and result.
        std::swap(slave_input, slave_output);
        std::swap(local_host_input, local_host_output);

        athread_dma_wait_value(&dma_reply, dma_reply_count);

        // Synchronize.
        if (i != args.iterations - 1)
            athread_ssync_array();
    }

    // Write the final result.
    slave_input.dma_put_content_to(local_host_input);
}
