#include "stencil/boundary_matrix.hpp"
#include "stencil_slave.hpp"

#include <algorithm>
#include <sys/cdefs.h>
#include <utility>

namespace {
struct RowTag { };
struct ColTag { };
constexpr RowTag row_tag;
constexpr ColTag col_tag;

template <class ElemTy, unsigned BoundaryWidth, unsigned BoundaryHeight>
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
        _data((ElemTy*) ldm_malloc(inner_data_mem_size())) {
        std::generate(std::begin(boundaries), std::begin(boundaries) + 2, [this] {
            return (ElemTy*) ldm_malloc(boundary_mem_size(BoundaryPosition::Top));
        });
        std::generate(std::begin(boundaries) + 2, std::end(boundaries), [this] {
            return (ElemTy*) ldm_malloc(boundary_mem_size(BoundaryPosition::Left));
        });
    }

    SlaveBlock(SlaveBlock&& other) noexcept :
        _rows(other._rows), _cols(other._cols), _data(other._data) {
        std::copy(
            std::begin(other.boundaries),
            std::end(other.boundaries),
            std::begin(boundaries)
        );
        std::fill(std::begin(other.boundaries), std::end(other.boundaries), nullptr);
        other._data = nullptr;
        other._rows = other._cols = 0;
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

            other._data = nullptr;
            std::fill(std::begin(other.boundaries), std::end(other.boundaries), nullptr);
            other._rows = other._cols = 0;
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

    template <int Col>
    auto data_at(unsigned row) const -> ElemTy& {
        if constexpr (Col >= 0) {
            if constexpr (Col < (int) BoundaryWidth) {
                return left_inner_boundary_at(row, Col);
            } else {
                return data_at(row, Col);
            }
        } else {
            if constexpr (Col >= -(int) BoundaryWidth) {
                return right_inner_boundary_at(row, Col + (int) BoundaryWidth);
            } else {
                return data_at(row, (int) _cols + Col);
            }
        }
    }

    auto data_at(unsigned row, unsigned col) const -> ElemTy& {
        return _data[row * (_cols - BoundaryWidth * 2) + col - BoundaryWidth];
    }

    template <int Row, int Col>
    auto data_with_boundary_at() const -> ElemTy& {
        if constexpr (Row >= 0) {
            if constexpr (Row < (int) BoundaryHeight) {
                // In top boundary, different Col have the same subscript calculation logic.
                return top_boundary_at(
                    Row,
                    Col < 0 ? (int) _cols + Col + (int) BoundaryWidth : Col - BoundaryWidth
                );
            } else {
                constexpr unsigned r = Row - BoundaryHeight;
                // We are accessing core element area.
                if constexpr (Col >= 0 && Col < (int) BoundaryWidth) {
                    return left_boundary_at(r, Col);
                } else if constexpr (Col < 0 && Col >= -(int) BoundaryWidth) {
                    return right_boundary_at(r, Col + (int) BoundaryWidth);
                } else {
                    // Fallback to data_at()
                    // clang-format off
                    return data_at<  //
                        /*Col=*/(Col < 0 ? Col + (int) BoundaryWidth : Col - (int) BoundaryWidth)
                    >(/*row=*/r);
                    // clang-format on
                }
            }
        } else {
            if constexpr (Row >= -(int) BoundaryHeight) {
                // In bottom boundary, different Col have the same subscript calculation logic.
                return bottom_boundary_at(
                    Row + (int) BoundaryHeight,
                    Col < 0 ? (int) _cols + Col + (int) BoundaryWidth : Col - BoundaryWidth
                );
            } else {
                // We are accessing core element area.
                if constexpr (Col >= 0 && Col < (int) BoundaryWidth) {
                    return left_boundary_at(Row + (int) BoundaryHeight + (int) _rows, Col);
                } else if constexpr (Col < 0 && Col >= -(int) BoundaryWidth) {
                    return right_boundary_at(
                        Row + (int) BoundaryHeight + (int) _rows,
                        Col + (int) BoundaryWidth
                    );
                } else {
                    // Fallback to data_at()
                    // clang-format off
                    return data_at<
                        /*Col=*/(Col < 0 ? Col + (int) BoundaryWidth : Col - (int) BoundaryWidth)
                    >(/*row=*/Row + (int) BoundaryHeight + (int) _rows);
                    // clang-format on
                }
            }
        }
    }

    template <int Row>
    auto data_with_boundary_at(RowTag, unsigned col) const -> ElemTy& {
        if constexpr (Row >= 0) {
            if constexpr (Row < (int) BoundaryHeight) {
                return top_boundary_at(Row, col - BoundaryWidth);
            } else {
                return data_at(Row - BoundaryHeight, col - BoundaryWidth);
            }
        } else {
            if constexpr (Row >= -(int) BoundaryHeight) {
                return bottom_boundary_at(Row + (int) BoundaryHeight, col - BoundaryWidth);
            } else {
                return data_at(Row + (int) BoundaryHeight + (int) _rows, col - BoundaryWidth);
            }
        }
    }

    template <int Col>
    auto data_with_boundary_at(unsigned row, ColTag) const -> ElemTy& {
        if constexpr (Col >= 0) {
            if constexpr (Col < (int) BoundaryWidth) {
                return left_boundary_at(row - BoundaryHeight, Col);
            } else {
                return data_at<Col - (int) BoundaryWidth>(row - BoundaryHeight);
            }
        } else {
            if constexpr (Col >= -(int) BoundaryWidth) {
                return right_boundary_at(row - BoundaryHeight, Col + (int) BoundaryWidth);
            } else {
                return data_at<Col + (int) BoundaryWidth>(row - BoundaryHeight);
            }
        }
    }

    auto top_boundary_at(unsigned row, unsigned col) const -> ElemTy& {
        return boundaries[Top][row * _cols + col];
    }

    auto bottom_boundary_at(unsigned row, unsigned col) const -> ElemTy& {
        return boundaries[Bottom][row * _cols + col];
    }

    auto left_boundary_at(unsigned row, unsigned col) const -> ElemTy& {
        return boundaries[Left][row * BoundaryWidth + col];
    }

    auto left_inner_boundary_at(unsigned row, unsigned col) const -> ElemTy& {
        return boundaries[LeftIn][row * BoundaryWidth + col];
    }

    auto right_boundary_at(unsigned row, unsigned col) const -> ElemTy& {
        return boundaries[Right][row * BoundaryWidth + col];
    }

    auto right_inner_boundary_at(unsigned row, unsigned col) const -> ElemTy& {
        return boundaries[RightIn][row * BoundaryWidth + col];
    }

    void dma_get_content_from(BoundaryMatrixView<ElemTy>& src_block) const {
        unsigned const vec_block_size = _cols - 2 * BoundaryWidth;

        // main body
        athread_dma_get_stride(
            &data_at(0, BoundaryWidth),
            &src_block.elem_at(0, BoundaryWidth),
            inner_data_mem_size(),
            vec_block_size * sizeof(ElemTy),
            (src_block.internal_data_stride() - vec_block_size) * sizeof(ElemTy)
        );

        // inner left
        athread_dma_get_stride(
            &left_inner_boundary_at(0, 0),
            &src_block.elem_at(0, 0),
            boundary_mem_size(BoundaryPosition::LeftIn),
            BoundaryWidth * sizeof(ElemTy),
            (src_block.internal_data_stride() - BoundaryWidth) * sizeof(ElemTy)
        );

        // inner right
        athread_dma_get_stride(
            &right_inner_boundary_at(0, 0),
            &src_block.elem_at(0, src_block.width() - BoundaryWidth),
            boundary_mem_size(BoundaryPosition::RightIn),
            BoundaryWidth * sizeof(ElemTy),
            (src_block.internal_data_stride() - BoundaryWidth) * sizeof(ElemTy)
        );
    }

    void dma_iget_top_boundary_from(  //
        BoundaryMatrixView<ElemTy>& src_block,
        athread_rply_t* reply
    ) const {
        if (BoundaryHeight == 1) {
            athread_dma_iget(
                &top_boundary_at(0, 0),
                &src_block.elem_with_boundary_at(0, BoundaryWidth),
                boundary_mem_size(BoundaryPosition::Top),
                reply
            );
        } else {
            athread_dma_iget_stride(
                &top_boundary_at(0, 0),
                &src_block.elem_with_boundary_at(0, BoundaryWidth),
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
        if (BoundaryHeight == 1) {
            athread_dma_iget(
                &bottom_boundary_at(0, 0),
                &src_block.elem_with_boundary_at(
                    src_block.height_with_boundary() - src_block.boundary_height(),
                    BoundaryWidth
                ),
                boundary_mem_size(BoundaryPosition::Bottom),
                reply
            );
        } else {
            athread_dma_iget_stride(
                &bottom_boundary_at(0, 0),
                &src_block.elem_with_boundary_at(
                    src_block.height_with_boundary() - src_block.boundary_height(),
                    BoundaryWidth
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
            BoundaryWidth * sizeof(ElemTy),
            (src_block.internal_data_stride() - BoundaryWidth) * sizeof(ElemTy),
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
            BoundaryWidth * sizeof(ElemTy),
            (src_block.internal_data_stride() - BoundaryWidth) * sizeof(ElemTy),
            reply
        );
    }

    void dma_iget_boundary_from(  //
        BoundaryMatrixView<ElemTy>& src_block,
        athread_rply_t* reply,
        unsigned& reply_count
    ) const {
        dma_iget_top_boundary_from(src_block, reply);
        dma_iget_bottom_boundary_from(src_block, reply);
        dma_iget_left_boundary_from(src_block, reply);
        dma_iget_right_boundary_from(src_block, reply);

        reply_count += 4;
    }

    void dma_put_content_to(BoundaryMatrixView<ElemTy>& dst_block) const {
        unsigned const vec_block_size = _cols - 2 * BoundaryWidth;

        // main body
        athread_dma_put_stride(
            &dst_block.elem_at(0, BoundaryWidth),
            &data_at(0, BoundaryWidth),
            inner_data_mem_size(),
            vec_block_size * sizeof(ElemTy),
            (dst_block.internal_data_stride() - vec_block_size) * sizeof(ElemTy)
        );

        // inner left
        athread_dma_put_stride(
            &dst_block.elem_at(0, 0),
            &left_inner_boundary_at(0, 0),
            boundary_mem_size(BoundaryPosition::LeftIn),
            BoundaryWidth * sizeof(ElemTy),
            (dst_block.internal_data_stride() - BoundaryWidth) * sizeof(ElemTy)
        );

        // inner right
        athread_dma_put_stride(
            &dst_block.elem_at(0, dst_block.width() - BoundaryWidth),
            &right_inner_boundary_at(0, 0),
            boundary_mem_size(BoundaryPosition::RightIn),
            BoundaryWidth * sizeof(ElemTy),
            (dst_block.internal_data_stride() - BoundaryWidth) * sizeof(ElemTy)
        );
    }

    void dma_iput_top_boundary_to(  //
        BoundaryMatrixView<ElemTy>& dst_block,
        athread_rply_t* reply
    ) const {
        if constexpr (BoundaryHeight == 1) {
            athread_dma_iput(
                &dst_block.elem_at(0, BoundaryWidth),
                &data_at(0, BoundaryWidth),
                (_cols - BoundaryWidth * 2) * sizeof(ElemTy),
                reply
            );
        } else {
            unsigned const vec_block_size = _cols - BoundaryWidth * 2;
            athread_dma_iput_stride(
                &dst_block.elem_at(0, BoundaryWidth),
                &data_at(0, BoundaryWidth),
                vec_block_size * BoundaryHeight * sizeof(ElemTy),
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
        if constexpr (BoundaryHeight == 1) {
            athread_dma_iput(
                &dst_block.elem_at(_rows - BoundaryHeight, BoundaryWidth),
                &data_at(_rows - BoundaryHeight, BoundaryWidth),
                (_cols - BoundaryWidth * 2) * sizeof(ElemTy),
                reply
            );
        } else {
            unsigned const vec_block_size = _cols - BoundaryWidth * 2;
            athread_dma_iput_stride(
                &dst_block.elem_at(_rows - BoundaryHeight, BoundaryWidth),
                &data_at(_rows - BoundaryHeight, BoundaryWidth),
                vec_block_size * BoundaryHeight * sizeof(ElemTy),
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
            BoundaryWidth * sizeof(ElemTy),
            (dst_block.internal_data_stride() - BoundaryWidth) * sizeof(ElemTy),
            reply
        );
    }

    void dma_iput_right_boundary_to(  //
        BoundaryMatrixView<ElemTy>& dst_block,
        athread_rply_t* reply
    ) const {
        athread_dma_iput_stride(
            &dst_block.elem_at(0, _cols - BoundaryWidth),
            &right_inner_boundary_at(0, 0),
            boundary_mem_size(BoundaryPosition::RightIn),
            BoundaryWidth * sizeof(ElemTy),
            (dst_block.internal_data_stride() - BoundaryWidth) * sizeof(ElemTy),
            reply
        );
    }

private:
    unsigned _rows;
    unsigned _cols;
    ElemTy* _data;
    ElemTy* boundaries[BoundaryCount];

    auto inner_data_mem_size() const -> std::size_t {
        return sizeof(ElemTy) * _rows * (_cols - BoundaryWidth * 2);
    }

    auto boundary_mem_size(BoundaryPosition position) const -> std::size_t {
        switch (position) {
        case Top:
        case Bottom:
            return _cols * BoundaryHeight * sizeof(ElemTy);
        case Left:
        case LeftIn:
        case Right:
        case RightIn:
            return _rows * BoundaryWidth * sizeof(ElemTy);
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

template <class ElemTy, unsigned BoundaryWidth, unsigned BoundaryHeight>
class StencilImpl {
public:
    explicit StencilImpl(Arguments& args) :
        args(args),
        local_host_input(args.input.block_subview(args.block_size, args.block_size, _ROW, _COL)),
        local_host_output(args.output.block_subview(args.block_size, args.block_size, _ROW, _COL)),
        slave_input(local_host_input),
        slave_output(local_host_output),
        width(local_host_input.width()),
        height(local_host_input.height()) { }

    void perform_iteration() {
        if (slave_input.empty() || slave_output.empty()) {
            return;
        }

        // Load content into LDM.
        slave_input.dma_get_content_from(local_host_input);

        athread_rply_t dma_reply = 0;
        unsigned dma_reply_count = 0;

        for (unsigned i = 0; i != args.iterations; ++i) {
            // Load boundaries.
            slave_input.dma_iget_boundary_from(local_host_input, &dma_reply, dma_reply_count);

            // Compute core elements.
            compute_core_elems(std::make_integer_sequence<unsigned, BoundaryWidth>());

            // Wait for boundaries.
            athread_dma_wait_value(&dma_reply, dma_reply_count);

            compute_left_top_corner(
                std::make_integer_sequence<unsigned, BoundaryWidth>(),
                std::make_integer_sequence<unsigned, BoundaryHeight>()
            );

            compute_left_boundary(std::make_integer_sequence<unsigned, BoundaryWidth>());

            compute_left_bottom_corner(
                std::make_integer_sequence<unsigned, BoundaryWidth>(),
                std::make_integer_sequence<unsigned, BoundaryHeight>()
            );

            // Write left boundary to the host memory.
            slave_output.dma_iput_left_boundary_to(local_host_output, &dma_reply);

            compute_right_top_corner(
                std::make_integer_sequence<unsigned, BoundaryWidth>(),
                std::make_integer_sequence<unsigned, BoundaryHeight>()
            );

            compute_right_boundary(std::make_integer_sequence<unsigned, BoundaryWidth>());

            compute_right_bottom_corner(
                std::make_integer_sequence<unsigned, BoundaryWidth>(),
                std::make_integer_sequence<unsigned, BoundaryHeight>()
            );

            // Write right boundary to the host memory.
            slave_output.dma_iput_right_boundary_to(local_host_output, &dma_reply);

            compute_top_boundary(
                std::make_integer_sequence<unsigned, BoundaryWidth>(),
                std::make_integer_sequence<unsigned, BoundaryHeight>()
            );

            // Write top boundary to the host memory.
            slave_output.dma_iput_top_boundary_to(local_host_output, &dma_reply);

            compute_bottom_boundary(
                std::make_integer_sequence<unsigned, BoundaryWidth>(),
                std::make_integer_sequence<unsigned, BoundaryHeight>()
            );

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

private:
    Arguments& args;
    BoundaryMatrixView<ElemTy> local_host_input, local_host_output;
    SlaveBlock<ElemTy, BoundaryWidth, BoundaryHeight> slave_input, slave_output;
    unsigned width, height;

    static constexpr ElemTy avg = ElemTy(1) / ElemTy((BoundaryWidth + BoundaryHeight) * 2);

    template <int Col, unsigned... BWIs, unsigned... BHIs>
    // clang-format off
    void core_stencil_computation(
        unsigned row,
        std::integer_sequence<unsigned, BWIs...>,
        std::integer_sequence<unsigned, BHIs...>
    ) {
        // clang-format on
        auto const sum =
            // Loop [Col - BoundaryWidth, Col): left neighbors.
            (... + (slave_input.template data_at</*Col=*/Col - BWIs - 1>(row)))
            // Loop [Col + 1, Col + BoundaryWidth + 1): right neighbors.
            + (... + (slave_input.template data_at</*Col=*/Col + 1 + BWIs>(row)))
            // Loop [row - BoundaryHeight, row): top neighbors.
            + (... + (slave_input.template data_at<Col>(/*row=*/row - BHIs - 1)))
            // Loop [row + 1, row + BoundaryHeight + 1): bottom neighbors.
            + (... + (slave_input.template data_at<Col>(/*row=*/row + BHIs + 1)));

        slave_output.template data_at<Col>(row) = sum * avg;
    }

    template <unsigned... BWIs, unsigned... BHIs>
    // clang-format off
    void core_stencil_computation(
        unsigned row,
        unsigned col,
        std::integer_sequence<unsigned, BWIs...>,
        std::integer_sequence<unsigned, BHIs...>
    ) {
        // clang-format on
        auto const sum =
            // Loop [Col - BoundaryWidth, Col): left neighbors.
            (... + (slave_input.data_at(row, col - BWIs - 1)))
            // Loop [Col + 1, Col + BoundaryWidth + 1): right neighbors.
            + (... + (slave_input.data_at(row, col + 1 + BWIs)))
            // Loop [row - BoundaryHeight, row): top neighbors.
            + (... + (slave_input.data_at(row - BHIs - 1, col)))
            // Loop [row + 1, row + BoundaryHeight + 1): bottom neighbors.
            + (... + (slave_input.data_at(row + BHIs + 1, col)));

        slave_output.data_at(row, col) = sum * avg;
    }

    template <unsigned... BWIs>
    void compute_core_elems(std::integer_sequence<unsigned, BWIs...>) {
        for (unsigned row = BoundaryHeight; row != height - BoundaryHeight; ++row) {
            // For core elements whose left neighbor may be stored in LeftIn buffer.
            // TODO: Note that if the neighborhood is too large, an element's left and right
            // neighbors may be saved in LeftIn and RightIn respectively in the same time. We don't
            // consider this case here.
            (core_stencil_computation</*Col=*/BWIs + BoundaryWidth>(
                 row,
                 std::make_integer_sequence<unsigned, BoundaryWidth>(),
                 std::make_integer_sequence<unsigned, BoundaryHeight>()
             ),
             ...);

            for (unsigned col = BoundaryWidth * 2; col != width - BoundaryWidth * 2; ++col) {
                core_stencil_computation(
                    row,
                    col,
                    std::make_integer_sequence<unsigned, BoundaryWidth>(),
                    std::make_integer_sequence<unsigned, BoundaryHeight>()
                );
            }

            (core_stencil_computation</*Col=*/(int) BWIs - (int) BoundaryWidth * 2>(
                 row,
                 std::make_integer_sequence<unsigned, BoundaryWidth>(),
                 std::make_integer_sequence<unsigned, BoundaryHeight>()
             ),
             ...);
        }
    }

    template <int Row, int Col, unsigned... BWIs, unsigned... BHIs>
    // clang-format off
    void corner_stencil_computation(
        std::integer_sequence<unsigned, BWIs...>,
        std::integer_sequence<unsigned, BHIs...>
    ) {
        // clang-format on
        auto const sum =
            // Loop [Col - BoundaryWidth, Col): left neighbors.
            (... + (slave_input.template data_with_boundary_at<Row, Col - BWIs - 1>()))
            // Loop [Col + 1, Col + BoundaryWidth + 1): right neighbors.
            + (... + (slave_input.template data_with_boundary_at<Row, Col + 1 + BWIs>()))
            // Loop [Row - BoundaryHeight, Row): top neighbors.
            + (... + (slave_input.template data_with_boundary_at<Row - BHIs - 1, Col>()))
            // Loop [Row + 1, Row + BoundaryHeight + 1): bottom neighbors.
            + (... + (slave_input.template data_with_boundary_at<Row + BHIs + 1, Col>()));

        slave_output.template data_with_boundary_at<Row, Col>() = sum * avg;
    }

    template <int Row, int... Cols>
    void compute_corner_impl() {
        (...,
         corner_stencil_computation<Row, Cols>(
             std::make_integer_sequence<unsigned, BoundaryWidth>(),
             std::make_integer_sequence<unsigned, BoundaryHeight>()
         ));
    }

    template <unsigned... BWIs, unsigned... BHIs>
    // clang-format off
    void compute_left_top_corner(
        std::integer_sequence<unsigned, BWIs...>,
        std::integer_sequence<unsigned, BHIs...>
    ) {
        // clang-format on
        (..., compute_corner_impl<BHIs + BoundaryHeight, (BWIs + BoundaryWidth)...>());
    }

    template <unsigned... BWIs, unsigned... BHIs>
    // clang-format off
    void compute_left_bottom_corner(
        std::integer_sequence<unsigned, BWIs...>,
        std::integer_sequence<unsigned, BHIs...>
    ) {
        // clang-format on
        (..., compute_corner_impl<BHIs - (int) BoundaryHeight * 2, (BWIs + BoundaryWidth)...>());
    }

    template <unsigned... BWIs, unsigned... BHIs>
    // clang-format off
    void compute_right_top_corner(
        std::integer_sequence<unsigned, BWIs...>,
        std::integer_sequence<unsigned, BHIs...>
    ) {
        // clang-format on
        (...,
         compute_corner_impl<BHIs + BoundaryHeight, ((int) BWIs - (int) BoundaryWidth * 2)...>());
    }

    template <unsigned... BWIs, unsigned... BHIs>
    // clang-format off
    void compute_right_bottom_corner(
        std::integer_sequence<unsigned, BWIs...>,
        std::integer_sequence<unsigned, BHIs...>
    ) {
        // clang-format on
        (...,
         compute_corner_impl<  //
             (int) BHIs - (int) BoundaryHeight * 2,
             ((int) BWIs - (int) BoundaryWidth * 2)...>());
    }

    template <int Col, unsigned... BWIs, unsigned... BHIs>
    // clang-format off
    void boundary_stencil_computation(
        unsigned row,
        ColTag,
        std::integer_sequence<unsigned, BWIs...>,
        std::integer_sequence<unsigned, BHIs...>
    ) {
        // clang-format on
        auto const sum =
            // Loop [Col - BoundaryWidth, Col): left neighbors.
            (... + (slave_input.template data_with_boundary_at<Col - BWIs - 1>(row, col_tag)))
            // Loop [Col + 1, Col + BoundaryWidth + 1): right neighbors.
            + (... + (slave_input.template data_with_boundary_at<Col + 1 + BWIs>(row, col_tag)))
            // Loop [Row - BoundaryHeight, Row): top neighbors.
            + (... + (slave_input.template data_with_boundary_at<Col>(row - BHIs - 1, col_tag)))
            // Loop [Row + 1, Row + BoundaryHeight + 1): bottom neighbors.
            + (... + (slave_input.template data_with_boundary_at<Col>(row + BHIs + 1, col_tag)));

        slave_output.template data_with_boundary_at<Col>(row, col_tag) = sum * avg;
    }

    template <int Row, unsigned... BWIs, unsigned... BHIs>
    // clang-format off
    void boundary_stencil_computation(
        RowTag,
        unsigned col,
        std::integer_sequence<unsigned, BWIs...>,
        std::integer_sequence<unsigned, BHIs...>
    ) {
        // clang-format on
        auto const sum =
            // Loop [Col - BoundaryWidth, Col): left neighbors.
            (... + (slave_input.template data_with_boundary_at<Row>(row_tag, col - BWIs - 1)))
            // Loop [Col + 1, Col + BoundaryWidth + 1): right neighbors.
            + (... + (slave_input.template data_with_boundary_at<Row>(row_tag, col + 1 + BWIs)))
            // Loop [Row - BoundaryHeight, Row): top neighbors.
            + (... + (slave_input.template data_with_boundary_at<Row - BHIs - 1>(row_tag, col)))
            // Loop [Row + 1, Row + BoundaryHeight + 1): bottom neighbors.
            + (... + (slave_input.template data_with_boundary_at<Row + BHIs + 1>(row_tag, col)));

        slave_output.template data_with_boundary_at<Row>(row_tag, col) = sum * avg;
    }

    template <unsigned... BWIs>
    void compute_left_boundary(std::integer_sequence<unsigned, BWIs...>) {
        for (unsigned row = BoundaryHeight * 2; row != height; ++row) {
            (...,
             boundary_stencil_computation</*Col=*/BWIs + BoundaryWidth>(
                 row,
                 col_tag,
                 std::make_integer_sequence<unsigned, BoundaryWidth>(),
                 std::make_integer_sequence<unsigned, BoundaryHeight>()
             ));
        }
    }

    template <unsigned... BWIs>
    void compute_right_boundary(std::integer_sequence<unsigned, BWIs...>) {
        for (unsigned row = BoundaryHeight * 2; row != height; ++row) {
            (...,
             boundary_stencil_computation</*Col=*/(int) BWIs - (int) BoundaryWidth * 2>(
                 row,
                 col_tag,
                 std::make_integer_sequence<unsigned, BoundaryWidth>(),
                 std::make_integer_sequence<unsigned, BoundaryHeight>()
             ));
        }
    }

    template <unsigned... BWIs, unsigned... BHIs>
    // clang-format off
    void compute_top_boundary(
        std::integer_sequence<unsigned, BWIs...>,
        std::integer_sequence<unsigned, BHIs...>
    ) {
        // clang-format on
        (..., compute_corner_impl<BHIs + BoundaryHeight, (BWIs + BoundaryWidth * 2)...>());

        for (unsigned col = BoundaryWidth * 3; col != width - BoundaryWidth; ++col) {
            (...,
             boundary_stencil_computation</*Row=*/BHIs + BoundaryHeight>(
                 row_tag,
                 col,
                 std::make_integer_sequence<unsigned, BoundaryWidth>(),
                 std::make_integer_sequence<unsigned, BoundaryHeight>()
             ));
        }

        (...,
         compute_corner_impl<  //
             BHIs + BoundaryHeight,
             ((int) BWIs - (int) BoundaryWidth * 3)...>());
    }

    template <unsigned... BWIs, unsigned... BHIs>
    // clang-format off
    void compute_bottom_boundary(
        std::integer_sequence<unsigned, BWIs...>,
        std::integer_sequence<unsigned, BHIs...>
    ) {
        // clang-format on
        (...,
         compute_corner_impl<  //
             (int) BHIs - (int) BoundaryHeight * 2,
             (BWIs + BoundaryWidth * 2)...>());

        for (unsigned col = BoundaryWidth * 3; col != width - BoundaryWidth; ++col) {
            (...,
             boundary_stencil_computation</*Row=*/(int) BHIs - (int) BoundaryHeight * 2>(
                 row_tag,
                 col,
                 std::make_integer_sequence<unsigned, BoundaryWidth>(),
                 std::make_integer_sequence<unsigned, BoundaryHeight>()
             ));
        }

        (...,
         compute_corner_impl<  //
             (int) BHIs - (int) BoundaryHeight * 2,
             ((int) BWIs - (int) BoundaryWidth * 3)...>());
    }
};

template <unsigned BoundaryWidth, unsigned BoundaryHeight>
void stencil_iterate_dma(Arguments& args) {
    StencilImpl<float, BoundaryWidth, BoundaryHeight> impl(args);
    impl.perform_iteration();
}
}  // namespace

void stencil_iterate_dma_static_unroll(Arguments* _args) {
    Arguments args;
    athread_dma_get(&args, _args, sizeof(Arguments));

    stencil_iterate_dma<1, 1>(args);
}
