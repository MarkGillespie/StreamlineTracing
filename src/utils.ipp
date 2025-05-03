template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& data) {
    out << "{";
    for (size_t i = 0; i < std::min(data.size(), (size_t)5); ++i) {
        out << data[i];
        if (i + 1 < data.size()) out << ", ";
    }
    if (data.size() >= 5) out << "... (length " << data.size() << ")";
    out << "}";
    return out;
}

template <typename T, size_t N>
std::ostream& operator<<(std::ostream& out, const std::array<T, N>& data) {
    out << "{";
    for (size_t i = 0; i < std::min(N, (size_t)5); ++i) {
        out << data[i];
        if (i + 1 < N) out << ", ";
    }
    if (N >= 5) out << "... (length " << N << ")";
    out << "}";
    return out;
}

template <typename T, typename S>
std::ostream& operator<<(std::ostream& out, const std::pair<T, S>& data) {
    out << "(" << std::get<0>(data) << ", " << std::get<1>(data) << ")";
    return out;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::set<T>& data) {
    out << "{";
    size_t iPrinted = 0;
    for (const T& element : data) {
        out << element;
        if (iPrinted + 1 < data.size()) out << ", ";
        if (iPrinted++ >= 4) break;
    }
    if (data.size() >= 5) out << "... (length " << data.size() << ")";
    out << "}";
    return out;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const Eigen::Triplet<T>& data) {
    out << "(row: " << data.row() << ", col: " << data.col()
        << ", value: " << data.value() << ")";
    return out;
}
