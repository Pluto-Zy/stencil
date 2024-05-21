#include "stencil/boundary_matrix.hpp"
#include "stencil_slave.hpp"

#include <algorithm>
#include <cstring>
#include <iterator>

namespace {
enum class RMAReply {
    Top,
    Bottom,
    Left,
    Right,
};

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
        return _data[row * _cols + col];
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
        // main body
        athread_dma_get_stride(
            &data_at(0, 0),
            &src_block.elem_at(0, 0),
            inner_data_mem_size(),
            _cols * sizeof(ElemTy),
            (src_block.internal_data_stride() - _cols) * sizeof(ElemTy)
        );

        for (unsigned r = 0; r != _rows; ++r) {
            // inner left
            std::memcpy(
                &left_inner_boundary_at(r, 0),
                &data_at(r, 0),
                _boundary_width * sizeof(ElemTy)
            );

            // inner right
            std::memcpy(
                &right_inner_boundary_at(r, 0),
                &data_at(r, _cols - _boundary_width),
                _boundary_width * sizeof(ElemTy)
            );
        }
    }

    void generate_boundary() const {
        if (_COL == 0) {
            std::fill_n(&left_boundary_at(0, 0), _rows, ElemTy(1));
        }

        if (_COL == 7) {
            std::fill_n(&right_boundary_at(0, 0), _rows, ElemTy(1));
        }

        if (_ROW == 0) {
            std::fill_n(&top_boundary_at(0, 0), _cols, ElemTy(0));
        }

        if (_ROW == 7) {
            std::fill_n(&bottom_boundary_at(0, 0), _cols, ElemTy(0));
        }
    }

    void dma_put_content_to(BoundaryMatrixView<ElemTy>& dst_block) const {
        for (unsigned r = 0; r != _rows; ++r) {
            // inner left
            std::memcpy(
                &data_at(r, 0),
                &left_inner_boundary_at(r, 0),
                _boundary_width * sizeof(ElemTy)
            );

            // inner right
            std::memcpy(
                &data_at(r, _cols - _boundary_width),
                &right_inner_boundary_at(r, 0),
                _boundary_width * sizeof(ElemTy)
            );
        }

        athread_dma_put_stride(
            &dst_block.elem_at(0, 0),
            &data_at(0, 0),
            inner_data_mem_size(),
            _cols * sizeof(ElemTy),
            (dst_block.internal_data_stride() - _cols) * sizeof(ElemTy)
        );
    }

    static auto slave_id_of(int row, int col) -> int {
        return row * 8 + col;
    }

    void rma_iput_top_boundary(athread_rply_t* send_reply, athread_rply_t* recv_reply) const {
        if (_ROW != 0) {
            athread_rma_iput(
                &data_at(0, 0),
                send_reply,
                boundary_mem_size(BoundaryPosition::Top),
                slave_id_of(_ROW - 1, _COL),
                &bottom_boundary_at(0, 0),
                recv_reply
            );
        }
    }

    void rma_iput_bottom_boundary(athread_rply_t* send_reply, athread_rply_t* recv_reply) const {
        if (_ROW != 7) {
            athread_rma_iput(
                &data_at(_rows - _boundary_height, 0),
                send_reply,
                boundary_mem_size(BoundaryPosition::Bottom),
                slave_id_of(_ROW + 1, _COL),
                &top_boundary_at(0, 0),
                recv_reply
            );
        }
    }

    void rma_iput_left_boundary(athread_rply_t* send_reply, athread_rply_t* recv_reply) const {
        if (_COL != 0) {
            athread_rma_iput(
                &left_inner_boundary_at(0, 0),
                send_reply,
                boundary_mem_size(BoundaryPosition::LeftIn),
                slave_id_of(_ROW, _COL - 1),
                &right_boundary_at(0, 0),
                recv_reply
            );
        }
    }

    void rma_iput_right_boundary(athread_rply_t* send_reply, athread_rply_t* recv_reply) const {
        if (_COL != 7) {
            athread_rma_iput(
                &right_inner_boundary_at(0, 0),
                send_reply,
                boundary_mem_size(BoundaryPosition::RightIn),
                slave_id_of(_ROW, _COL + 1),
                &left_boundary_at(0, 0),
                recv_reply
            );
        }
    }

    void rma_iput_boundaries(athread_rply_t* send_reply, athread_rply_t* recv_replies) const {
        rma_iput_bottom_boundary(send_reply, recv_replies + RMAReply::Top);
        rma_iput_top_boundary(send_reply, recv_replies + RMAReply::Bottom);
        rma_iput_right_boundary(send_reply, recv_replies + RMAReply::Left);
        rma_iput_left_boundary(send_reply, recv_replies + RMAReply::Right);
    }

private:
    unsigned _rows;
    unsigned _cols;
    unsigned _boundary_height;
    unsigned _boundary_width;
    ElemTy* _data;
    ElemTy* boundaries[BoundaryCount];

    auto inner_data_mem_size() const -> std::size_t {
        return sizeof(ElemTy) * _rows * _cols;
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

namespace neighbor1_impl {
void stencil_iterate_rma(Arguments& args) {
    // Divide the host matrix into blocks.
    auto local_host_input = args.input.block_subview(args.block_size, args.block_size, _ROW, _COL);
    auto local_host_output =
        args.output.block_subview(args.block_size, args.block_size, _ROW, _COL);

    // Allocate LDM.
    SlaveBlock slave_input(local_host_input);
    SlaveBlock slave_output(local_host_output);

    if (slave_input.empty() || slave_output.empty()) {
        return;
    }

    // Load content into LDM.
    slave_input.dma_get_content_from(local_host_input);
    slave_input.generate_boundary();
    slave_output.generate_boundary();

    athread_rply_t replies[10] {};

    // Iterations.
    for (unsigned i = 0; i != args.iterations; ++i) {
        athread_rply_t* send_reply = &replies[i % 2 * 5];
        athread_rply_t* recv_replies = send_reply + 1;

        // Send boundaries to other CPE.
        slave_input.rma_iput_boundaries(send_reply, recv_replies);

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

        // Wait for the top boundary.
        if (_ROW != 0) {
            athread_rma_wait_value(recv_replies + RMAReply::Top, 1);
            recv_replies[RMAReply::Top] = 0;
        }

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

        // Wait for the bottom bounary.
        if (_ROW != 7) {
            athread_rma_wait_value(recv_replies + RMAReply::Bottom, 1);
            recv_replies[RMAReply::Bottom] = 0;
        }

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

        // Wait for the left boundary.
        if (_COL != 0) {
            athread_rma_wait_value(recv_replies + RMAReply::Left, 1);
            recv_replies[RMAReply::Left] = 0;
        }

        {
            constexpr unsigned col = 0;
            unsigned row = 0;
            // Left-top corner.
            slave_output.left_inner_boundary_at(row, 0) = slave_output.data_at(row, col) = 0.25f
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
            slave_output.left_inner_boundary_at(row, 0) = slave_output.data_at(row, col) = 0.25f
                * (slave_input.left_inner_boundary_at(row - 1, 0)
                   + slave_input.left_boundary_at(row, 0) + slave_input.data_at(row, col + 1)
                   + slave_input.bottom_boundary_at(0, col));
        }

        // Wait for the right boundary.
        if (_COL != 7) {
            athread_rma_wait_value(recv_replies + RMAReply::Right, 1);
            recv_replies[RMAReply::Right] = 0;
        }

        {
            unsigned const col = slave_input.cols() - 1;
            unsigned row = 0;

            // Right-top corner.
            slave_output.right_inner_boundary_at(row, 0) = slave_output.data_at(row, col) = 0.25f
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
            slave_output.right_inner_boundary_at(row, 0) = slave_output.data_at(row, col) = 0.25f
                * (slave_input.right_inner_boundary_at(row - 1, 0)
                   + slave_input.data_at(row, col - 1) + slave_input.right_boundary_at(row, 0)
                   + slave_input.bottom_boundary_at(0, col));
        }

        // Swap the input and result.
        std::swap(slave_input, slave_output);
        std::swap(local_host_input, local_host_output);
    }

    // Write the final result.
    slave_input.dma_put_content_to(local_host_input);
}
}  // namespace neighbor1_impl
}  // namespace

void stencil_iterate_rma(Arguments* _args) {
    Arguments args;
    athread_dma_get(&args, _args, sizeof(Arguments));

    neighbor1_impl::stencil_iterate_rma(args);
}
