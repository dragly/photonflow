#ifndef RECTANGLE_H
#define RECTANGLE_H

#include "../core/geometry.h"

namespace photonflow {

class Rectangle
{
public:
    Rectangle();
    Rectangle(double x, double y, double width, double height);

    double x() const;
    double y() const;
    double width() const;
    double height() const;

private:
    double m_x = 0.0;
    double m_y = 0.0;
    double m_width = 0.0;
    double m_height = 0.0;
};

} // namespace

#endif // RECTANGLE_H
