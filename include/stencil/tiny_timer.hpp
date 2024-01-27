#ifndef TINY_TIMER_HPP
#define TINY_TIMER_HPP

#include <chrono>
#include <ostream>
#include <ratio>

namespace detail {
template <class Period>
struct DurationWrapper {
    std::chrono::duration<double, Period> duration;

    static constexpr auto unit(std::ratio<1>) -> const char* {
        return "s";
    }

    static constexpr auto unit(std::milli) -> const char* {
        return "ms";
    }

    static constexpr auto unit(std::micro) -> const char* {
        return "us";
    }

    static constexpr auto unit(std::nano) -> const char* {
        return "ns";
    }

    static constexpr auto unit() -> const char* {
        return unit(Period());
    }
};

template <class Period>
auto operator<<(std::ostream& out, DurationWrapper<Period> duration) -> std::ostream& {
    return out << duration.duration.count() << duration.unit();
}
}  // namespace detail

class TinyTimer {
    using Duration = std::chrono::steady_clock::duration;
    using TimePoint = std::chrono::steady_clock::time_point;

public:
    constexpr TinyTimer() : _start_time(), _total_duration() { }

    void start() {
        _start_time = std::chrono::steady_clock::now();
    }

    void pause() {
        auto const now = std::chrono::steady_clock::now();
        _total_duration += now - _start_time;
    }

    void restart() {
        _total_duration = {};
        start();
    }

    auto duration() const -> Duration {
        return _total_duration;
    }

    template <class Rep, class Period>
    auto duration() const -> std::chrono::duration<Rep, Period> {
        return std::chrono::duration_cast<std::chrono::duration<Rep, Period>>(_total_duration);
    }

    template <class Period>
    auto as() const -> detail::DurationWrapper<Period> {
        return { (std::chrono::duration<double, Period>) _total_duration };
    }

private:
    TimePoint _start_time;
    Duration _total_duration;
};

#endif  // TINY_TIMER_HPP
