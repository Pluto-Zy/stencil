#ifndef BOUNDARY_MATRIX
#define BOUNDARY_MATRIX

#include <cassert>
#include <cstdio>
#include <memory>
#include <tuple>
#include <type_traits>

enum class BoundaryPosition {
    Left,
    Right,
    Top,
    Bottom,
};

namespace detail {
/// Represents a matrix with additional boundaries for stencil computation.
///
/// x x x x x x
/// x x 0 0 x x
/// x x 0 0 x x
/// x x 0 0 x x
/// x x x x x x
///
/// Note that in this example, the matrix should be constructed with
///     - width: 2
///     - height: 3
///     - neighbor_width: 2
///     - neighbor_height: 1
template <class ElemTy, bool Owner>
class BoundaryMatrix {
public:
    using value_type = ElemTy;

    BoundaryMatrix() : _neighbor_width(0) { }

    /// Initialize the matrix with the given `width`, `height`, `neighbor_width` and
    /// `neighbor_height`. Note that `width` and `height` should be the size of the matrix
    /// *excluding* the width or height of the neighborhood.
    template <bool IsOwner = Owner, std::enable_if_t<IsOwner, int> = 0>
    BoundaryMatrix(
        std::size_t width,
        std::size_t height,
        unsigned neighbor_width,
        unsigned neighbor_height
    ) :
        _actual_width(width + neighbor_width * 2),
        _actual_height(height + neighbor_height * 2),
        _neighbor_width(neighbor_width),
        _neighbor_height(neighbor_height),
        _data_stride(_actual_width),
        _data(std::make_unique<ElemTy[]>(_actual_width * _actual_height)) { }

    template <bool IsOwner = Owner, std::enable_if_t<!IsOwner, int> = 0>
    BoundaryMatrix(
        std::size_t width,
        std::size_t height,
        unsigned neighbor_width,
        unsigned neighbor_height,
        std::size_t stride,
        ElemTy* data
    ) :
        _actual_width(width + neighbor_width * 2),
        _actual_height(height + neighbor_height * 2),
        _neighbor_width(neighbor_width),
        _neighbor_height(neighbor_height),
        _data_stride(stride),
        _data(data) { }

    auto width() const -> std::size_t {
        return _actual_width - _neighbor_width * 2;
    }

    auto height() const -> std::size_t {
        return _actual_height - _neighbor_height * 2;
    }

    auto width_with_boundary() const -> std::size_t {
        return _actual_width;
    }

    auto height_with_boundary() const -> std::size_t {
        return _actual_height;
    }

    auto boundary_width() const -> unsigned {
        return _neighbor_width;
    }

    auto boundary_height() const -> unsigned {
        return _neighbor_height;
    }

    auto internal_data_stride() const -> std::size_t {
        return _data_stride;
    }

    /// Returns the number of elements (excluding boundaries).
    auto elem_count() const -> std::size_t {
        return width() * height();
    }

    /// Returns the number of elements (including boundaries).
    auto elem_with_boundary_count() const -> std::size_t {
        return width_with_boundary() * height_with_boundary();
    }

    auto elem_with_boundary_at(std::size_t row, std::size_t col) const -> ElemTy& {
        assert(row < height_with_boundary() && col < width_with_boundary());
        return _data[row * _data_stride + col];
    }

    auto elem_at(std::size_t row, std::size_t col) const -> ElemTy& {
        assert(row < height() && col < width());
        return elem_with_boundary_at(row + _neighbor_height, col + _neighbor_width);
    }

    auto empty() const -> bool {
        return width() == 0 || height() == 0;
    }

    void fill_boundary(BoundaryPosition position, ElemTy value) {
        auto const [rb, re, cb, ce] = [&] {
            switch (position) {
            case BoundaryPosition::Left:
                return std::make_tuple(
                    std::size_t(0),
                    height_with_boundary(),
                    std::size_t(0),
                    std::size_t(boundary_width())
                );
            case BoundaryPosition::Right:
                return std::make_tuple(
                    std::size_t(0),
                    height_with_boundary(),
                    width_with_boundary() - boundary_width(),
                    width_with_boundary()
                );
            case BoundaryPosition::Top:
                return std::make_tuple(
                    std::size_t(0),
                    std::size_t(boundary_height()),
                    std::size_t(0),
                    width_with_boundary()
                );
            case BoundaryPosition::Bottom:
                return std::make_tuple(
                    height_with_boundary() - boundary_height(),
                    height_with_boundary(),
                    std::size_t(0),
                    width_with_boundary()
                );
            default:
                __builtin_unreachable();
            }
        }();

        for (std::size_t r = rb; r != re; ++r) {
            for (std::size_t c = cb; c != ce; ++c) {
                elem_with_boundary_at(r, c) = value;
            }
        }
    }

    auto borrow() const -> BoundaryMatrix<ElemTy, false> {
        // clang-format off
        return {
            width(),
            height(),
            boundary_width(),
            boundary_height(),
            internal_data_stride(),
            raw_data(),
        };
        // clang-format on
    }

    /// Create a subview of the matrix. It is like dividing the entire matrix into multiple blocks
    /// with a height of `block_height` and a width of `block_width`, and select the block at
    /// `row_idx` row and `col_idx` column. Note that the new block returned is a subview of the
    /// original matrix, and the underlying data is not copied. The neighborhood size of the new
    /// view is the same as the size of the original matrix.
    [[nodiscard]] auto block_subview(
        std::size_t block_width,
        std::size_t block_height,
        unsigned row_idx,
        unsigned col_idx
    ) const -> BoundaryMatrix<ElemTy, false> {
        std::size_t const subview_row_beg = row_idx * block_height;
        std::size_t const subview_col_beg = col_idx * block_width;

        if (subview_row_beg >= height() || subview_col_beg >= width()) {
            return {
                /*width=*/0,
                /*height=*/0,
                /*neighbor_width=*/boundary_width(),
                /*neightbor_height=*/boundary_height(),
                /*stride=*/internal_data_stride(),
                raw_data(),
            };
        } else {
            return {
                /*width=*/(std::min)(block_width, width() - subview_col_beg),
                /*height=*/(std::min)(block_height, height() - subview_row_beg),
                /*neighbor_width=*/boundary_width(),
                /*neighbor_height=*/boundary_height(),
                /*stride=*/internal_data_stride(),
                /*data=*/&elem_with_boundary_at(subview_row_beg, subview_col_beg),
            };
        }
    }

protected:
    std::size_t _actual_width;
    std::size_t _actual_height;
    unsigned _neighbor_width;
    unsigned _neighbor_height;
    std::size_t _data_stride;
    std::conditional_t<Owner, std::unique_ptr<ElemTy[]>, ElemTy*> _data;

    template <class E, bool O>
    friend class BoundaryMatrix;

    auto raw_data() const -> ElemTy* {
        if constexpr (Owner) {
            return _data.get();
        } else {
            return _data;
        }
    }
};
}  // namespace detail

template <class ElemTy>
using BoundaryMatrix = detail::BoundaryMatrix<ElemTy, true>;

template <class ElemTy>
using BoundaryMatrixView = detail::BoundaryMatrix<ElemTy, false>;

#endif  // BOUNDARY_MATRIX
