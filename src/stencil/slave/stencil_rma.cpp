#include "stencil_slave.hpp"
#include <algorithm>
#include <iterator>

namespace {
class HostBlock {
public:
    HostBlock(double* data, unsigned rows, unsigned cols) :
        _data(data), _rows(rows), _cols(cols), _stride(cols + 2), boundary(0) { }

    auto rows() const -> unsigned {
        return _rows;
    }

    auto cols() const -> unsigned {
        return _cols;
    }

    auto stride() const -> unsigned {
        return _stride;
    }

    auto boundary_width() const -> unsigned {
        return boundary;
    }

    auto empty() const -> bool {
        return _rows == 0 || _cols == 0;
    }

    auto data_at(unsigned row, unsigned col) const -> double& {
        return _data[row * _stride + col];
    }

    [[nodiscard]] auto divide_row_col(  //
        unsigned block_size,
        unsigned row_idx,
        unsigned col_idx
    ) const -> HostBlock {
        unsigned const block_row_beg = row_idx * block_size;
        unsigned const block_col_beg = col_idx * block_size;

        if (block_row_beg >= _rows || block_col_beg >= _cols) {
            return HostBlock(_data, 0, 0);
        } else {
            return HostBlock(
                _data + block_row_beg * _stride + block_col_beg,
                (std::min)(block_size, _rows - block_row_beg),
                (std::min)(block_size, _cols - block_col_beg),
                _stride
            );
        }
    }

    void expand_boundary(unsigned boundary_width) {
        _data = _data - boundary_width * _stride - boundary_width;
        _rows += 2 * boundary_width;
        _cols += 2 * boundary_width;
        boundary = boundary_width;
    }

private:
    double* _data;
    unsigned _rows;
    unsigned _cols;
    unsigned _stride;
    unsigned boundary;

    HostBlock(double* data, unsigned rows, unsigned cols, unsigned stride) :
        _data(data), _rows(rows), _cols(cols), _stride(stride), boundary(0) { }
};

class SlaveBlock {
    enum {
        Top,
        Bottom,
        LeftIn,
        Left,
        RightIn,
        Right,
        BoundaryCount,
    };

public:
    SlaveBlock(unsigned rows, unsigned cols) :
        _data((double*) ldm_malloc(sizeof(double) * rows * cols)), _rows(rows), _cols(cols) {
        std::generate(std::begin(boundaries), std::begin(boundaries) + 2, [this] {
            return (double*) ldm_malloc(sizeof(double) * _cols);
        });
        std::generate(std::begin(boundaries) + 2, std::end(boundaries), [this] {
            return (double*) ldm_malloc(sizeof(double) * _rows);
        });
    }

    SlaveBlock(SlaveBlock&& other) noexcept :
        _data(other._data), _rows(other._rows), _cols(other._cols) {
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
            || std::any_of(std::begin(boundaries), std::end(boundaries), [](double* ptr) {
                   return ptr == nullptr;
               });
    }

    auto data_at(unsigned row, unsigned col) const -> double& {
        return _data[row * _cols + col];
    }

    auto top_boundary_at(unsigned col) const -> double& {
        return boundaries[Top][col];
    }

    auto bottom_boundary_at(unsigned col) const -> double& {
        return boundaries[Bottom][col];
    }

    auto left_boundary_at(unsigned row) const -> double& {
        return boundaries[Left][row];
    }

    auto left_inner_boundary_at(unsigned row) const -> double& {
        return boundaries[LeftIn][row];
    }

    auto right_boundary_at(unsigned row) const -> double& {
        return boundaries[Right][row];
    }

    auto right_inner_boundary_at(unsigned row) const -> double& {
        return boundaries[RightIn][row];
    }

    void dma_get_content_from(const HostBlock& src_block) const {
        athread_dma_get_stride(
            &data_at(0, 0),
            &src_block.data_at(src_block.boundary_width(), src_block.boundary_width()),
            _rows * _cols * sizeof(double),
            _cols * sizeof(double),
            (src_block.stride() - _cols) * sizeof(double)
        );

        athread_dma_get_stride(
            &left_inner_boundary_at(0),
            &src_block.data_at(src_block.boundary_width(), src_block.boundary_width()),
            _rows * sizeof(double),
            sizeof(double),
            (src_block.stride() - 1) * sizeof(double)
        );

        athread_dma_get_stride(
            &right_inner_boundary_at(0),
            &src_block.data_at(
                src_block.boundary_width(),
                src_block.cols() - src_block.boundary_width() - 1
            ),
            _rows * sizeof(double),
            sizeof(double),
            (src_block.stride() - 1) * sizeof(double)
        );

        data_at(0, 0) = left_inner_boundary_at(0);
        data_at(0, _cols - 1) = right_inner_boundary_at(0);
        data_at(_rows - 1, 0) = left_inner_boundary_at(_rows - 1);
        data_at(_rows - 1, _cols - 1) = right_inner_boundary_at(_rows - 1);
    }

    void generate_boundary() const {
        if (_COL == 0) {
            std::fill_n(&left_boundary_at(0), _rows, 1.0);
        }

        if (_COL == 7) {
            std::fill_n(&right_boundary_at(0), _rows, 1.0);
        }

        if (_ROW == 0) {
            std::fill_n(&top_boundary_at(0), _cols, 0.0);
        }

        if (_ROW == 7) {
            std::fill_n(&bottom_boundary_at(0), _cols, 0.0);
        }
    }

    void dma_put_content_to(const HostBlock& dst_block) const {
        athread_dma_put_stride(
            &dst_block.data_at(dst_block.boundary_width(), dst_block.boundary_width()),
            &data_at(0, 0),
            _rows * _cols * sizeof(double),
            _cols * sizeof(double),
            (dst_block.stride() - _cols) * sizeof(double)
        );

        athread_dma_put_stride(
            &dst_block.data_at(dst_block.boundary_width(), dst_block.boundary_width()),
            &left_inner_boundary_at(0),
            _rows * sizeof(double),
            sizeof(double),
            (dst_block.stride() - 1) * sizeof(double)
        );

        athread_dma_put_stride(
            &dst_block.data_at(
                dst_block.boundary_width(),
                dst_block.cols() - dst_block.boundary_width() - 1
            ),
            &right_inner_boundary_at(0),
            _rows * sizeof(double),
            sizeof(double),
            (dst_block.stride() - 1) * sizeof(double)
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
                _cols * sizeof(double),
                slave_id_of(_ROW - 1, _COL),
                &bottom_boundary_at(0),
                recv_reply
            );
        }
    }

    void rma_iput_bottom_boundary(athread_rply_t* send_reply, athread_rply_t* recv_reply) const {
        if (_ROW != 7) {
            athread_rma_iput(
                &data_at(_rows - 1, 0),
                send_reply,
                _cols * sizeof(double),
                slave_id_of(_ROW + 1, _COL),
                &top_boundary_at(0),
                recv_reply
            );
        }
    }

    void rma_iput_left_boundary(athread_rply_t* send_reply, athread_rply_t* recv_reply) const {
        if (_COL != 0) {
            athread_rma_iput(
                &left_inner_boundary_at(0),
                send_reply,
                _rows * sizeof(double),
                slave_id_of(_ROW, _COL - 1),
                &right_boundary_at(0),
                recv_reply
            );
        }
    }

    void rma_iput_right_boundary(athread_rply_t* send_reply, athread_rply_t* recv_reply) const {
        if (_COL != 7) {
            athread_rma_iput(
                &right_inner_boundary_at(0),
                send_reply,
                _rows * sizeof(double),
                slave_id_of(_ROW, _COL + 1),
                &left_boundary_at(0),
                recv_reply
            );
        }
    }

    void rma_iput_boundaries(athread_rply_t* send_reply, athread_rply_t* recv_reply) const {
        rma_iput_top_boundary(send_reply, send_reply);
        rma_iput_left_boundary(send_reply, send_reply);
        rma_iput_bottom_boundary(send_reply, send_reply);
        rma_iput_right_boundary(send_reply, send_reply);
        (void) recv_reply;
    }

private:
    double* _data;
    double* boundaries[BoundaryCount];
    unsigned _rows;
    unsigned _cols;

    void free_resource() noexcept {
        if (_data) {
            ldm_free(_data, sizeof(double) * _rows * _cols);
        }

        std::for_each(std::begin(boundaries), std::begin(boundaries) + 2, [this](double* ptr) {
            if (ptr) {
                ldm_free(ptr, sizeof(double) * _cols);
            }
        });

        std::for_each(std::begin(boundaries) + 2, std::end(boundaries), [this](double* ptr) {
            if (ptr) {
                ldm_free(ptr, sizeof(double) * _rows);
            }
        });
    }
};
}  // namespace

void stencil_iterate_rma(Arguments* _args) {
    Arguments args;
    athread_dma_get(&args, _args, sizeof(Arguments));

    // Allocate LDM.
    SlaveBlock slave_input(args.block_size, args.block_size);
    SlaveBlock slave_output(args.block_size, args.block_size);

    if (slave_input.empty() || slave_output.empty()) {
        return;
    }

    // Create block for the whole matrix.
    HostBlock const host_input(args.input, args.matrix_size, args.matrix_size);
    HostBlock const host_output(args.output, args.matrix_size, args.matrix_size);
    HostBlock local_host_input = host_input.divide_row_col(args.block_size, _ROW, _COL);
    HostBlock local_host_output = host_output.divide_row_col(args.block_size, _ROW, _COL);
    // Expand the boundary.
    local_host_input.expand_boundary(1);
    local_host_output.expand_boundary(1);

    // Load content into LDM.
    slave_input.dma_get_content_from(local_host_input);
    slave_input.generate_boundary();
    slave_output.generate_boundary();

    unsigned const send_count = static_cast<unsigned>(_ROW != 0) + static_cast<unsigned>(_ROW != 7)
        + static_cast<unsigned>(_COL != 0) + static_cast<unsigned>(_COL != 7);
    // unsigned const recv_count = send_count;

    // Iterations.
    for (unsigned i = 0; i != args.iterations; ++i) {
        athread_rply_t send_reply = 0;
        athread_rply_t recv_reply = 0;

        // Send boundaries to other CPE.
        slave_input.rma_iput_boundaries(&send_reply, &recv_reply);

        // Compute contents.
        for (unsigned row = 1; row != slave_input.rows() - 1; ++row) {
            unsigned col = 1;
            slave_output.data_at(row, col) = 0.25
                * (slave_input.data_at(row - 1, col) + slave_input.left_inner_boundary_at(row)
                   + slave_input.data_at(row, col + 1) + slave_input.data_at(row + 1, col));

            for (++col; col != slave_input.cols() - 2; ++col) {
                slave_output.data_at(row, col) = 0.25
                    * (slave_input.data_at(row - 1, col) + slave_input.data_at(row, col - 1)
                       + slave_input.data_at(row, col + 1) + slave_input.data_at(row + 1, col));
            }

            slave_output.data_at(row, col) = 0.25
                * (slave_input.data_at(row - 1, col) + slave_input.data_at(row, col - 1)
                   + slave_input.right_inner_boundary_at(row) + slave_input.data_at(row + 1, col));
        }

        // Wait for boundaries.
        // athread_rma_wait_value(&recv_reply, recv_count);
        athread_rma_wait_value(&send_reply, send_count);

        {
            constexpr unsigned col = 0;
            unsigned row = 0;
            // Left-top corner.
            slave_output.left_inner_boundary_at(row) = slave_output.data_at(row, col) = 0.25
                * (slave_input.top_boundary_at(col) + slave_input.left_boundary_at(row)
                   + slave_input.data_at(row, col + 1)
                   + slave_input.left_inner_boundary_at(row + 1));

            // Left boundary.
            for (++row; row != slave_input.rows() - 1; ++row) {
                slave_output.left_inner_boundary_at(row) = 0.25
                    * (slave_input.left_inner_boundary_at(row - 1)
                       + slave_input.left_boundary_at(row) + slave_input.data_at(row, col + 1)
                       + slave_input.left_inner_boundary_at(row + 1));
            }

            // Left-bottom corner.
            slave_output.left_inner_boundary_at(row) = slave_output.data_at(row, col) = 0.25
                * (slave_input.left_inner_boundary_at(row - 1) + slave_input.left_boundary_at(row)
                   + slave_input.data_at(row, col + 1) + slave_input.bottom_boundary_at(col));
        }

        {
            unsigned const col = slave_input.cols() - 1;
            unsigned row = 0;

            // Right-top corner.
            slave_output.right_inner_boundary_at(row) = slave_output.data_at(row, col) = 0.25
                * (slave_input.top_boundary_at(col) + slave_input.data_at(row, col - 1)
                   + slave_input.right_boundary_at(row)
                   + slave_input.right_inner_boundary_at(row + 1));

            // Right boundary.
            for (++row; row != slave_input.rows() - 1; ++row) {
                slave_output.right_inner_boundary_at(row) = 0.25
                    * (slave_input.right_inner_boundary_at(row - 1)
                       + slave_input.data_at(row, col - 1) + slave_input.right_boundary_at(row)
                       + slave_input.right_inner_boundary_at(row + 1));
            }

            // Right-bottom corner.
            slave_output.right_inner_boundary_at(row) = slave_output.data_at(row, col) = 0.25
                * (slave_input.right_inner_boundary_at(row - 1) + slave_input.data_at(row, col - 1)
                   + slave_input.right_boundary_at(row) + slave_input.bottom_boundary_at(col));
        }

        // Top boundary.
        {
            unsigned col = 1;
            constexpr unsigned row = 0;

            slave_output.data_at(row, col) = 0.25
                * (slave_input.top_boundary_at(col) + slave_input.left_inner_boundary_at(row)
                   + slave_input.data_at(row, col + 1) + slave_input.data_at(row + 1, col));

            for (++col; col != slave_input.cols() - 2; ++col) {
                slave_output.data_at(row, col) = 0.25
                    * (slave_input.top_boundary_at(col) + slave_input.data_at(row, col - 1)
                       + slave_input.data_at(row, col + 1) + slave_input.data_at(row + 1, col));
            }

            slave_output.data_at(row, col) = 0.25
                * (slave_input.top_boundary_at(col) + slave_input.data_at(row, col - 1)
                   + slave_input.right_inner_boundary_at(row) + slave_input.data_at(row + 1, col));
        }

        // Bottom boundary.
        {
            unsigned const row = slave_input.rows() - 1;
            unsigned col = 1;

            slave_output.data_at(row, col) = 0.25
                * (slave_input.data_at(row - 1, col) + slave_input.left_inner_boundary_at(row)
                   + slave_input.data_at(row, col + 1) + slave_input.bottom_boundary_at(col));

            for (++col; col != slave_input.cols() - 2; ++col) {
                slave_output.data_at(row, col) = 0.25
                    * (slave_input.data_at(row - 1, col) + slave_input.data_at(row, col - 1)
                       + slave_input.data_at(row, col + 1) + slave_input.bottom_boundary_at(col));
            }

            slave_output.data_at(row, col) = 0.25
                * (slave_input.data_at(row - 1, col) + slave_input.data_at(row, col - 1)
                   + slave_input.right_inner_boundary_at(row)
                   + slave_input.bottom_boundary_at(col));
        }

        // Swap the input and result.
        std::swap(slave_input, slave_output);
        std::swap(local_host_input, local_host_output);

        if (_PEN == 0) {
            // printf("iteration %u\n", i);
        }

        athread_ssync_array();
    }

    // Write the final result.
    slave_input.dma_put_content_to(local_host_input);
}
