#ifndef OUTPUTCAPTURE_H
#define OUTPUTCAPTURE_H

#include <iostream>
#include <sstream>

namespace amr_wind_tests {

/** Class to handle output capture
 */
class CaptureOutput
{
public:
    CaptureOutput() { redirect_io(); }

    CaptureOutput(bool stdout, bool stderr) : m_stdout(stdout), m_stderr(stderr)
    {
        redirect_io();
    }

    ~CaptureOutput()
    {
        if (m_stdout && (m_cout_buf != nullptr)) {
            std::cout.rdbuf(m_cout_buf);
        }
        if (m_stderr && (m_cerr_buf != nullptr)) {
            std::cerr.rdbuf(m_cerr_buf);
        }
    }

    bool capture_stdout() const { return m_stdout; }
    bool capture_stderr() const { return m_stderr; }

    std::stringstream& stdout() { return m_cout; }
    std::stringstream& stderr() { return m_cerr; }

protected:
    void redirect_io()
    {
        if (m_stdout) {
            m_cout_buf = std::cout.rdbuf();
            std::cout.rdbuf(m_cout.rdbuf());
        }
        if (m_stderr) {
            m_cerr_buf = std::cerr.rdbuf();
            std::cerr.rdbuf(m_cerr.rdbuf());
        }
    }

    std::stringstream m_cout;
    std::stringstream m_cerr;

    std::streambuf* m_cout_buf{nullptr};
    std::streambuf* m_cerr_buf{nullptr};
    bool m_stdout{true};
    bool m_stderr{false};
};

} // namespace amr_wind_tests

#endif /* OUTPUTCAPTURE_H */
