#include "rectangle.h"

namespace photonflow {

Rectangle::Rectangle()
{

}

Rectangle::Rectangle(double x, double y, double width, double height)
    : m_x(x)
    , m_y(y)
    , m_width(width)
    , m_height(height)
{

}

double Rectangle::x() const
{
    return m_x;
}

double Rectangle::y() const
{
    return m_y;
}

double Rectangle::width() const
{
    return m_width;
}

double Rectangle::height() const
{
    return m_height;
}
} // namespace
