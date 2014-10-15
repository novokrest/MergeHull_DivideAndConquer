#include "stdafx.h"

#include "convex_hull.h"

using geom::structures::point_type;
using geom::structures::contour_type;
using geom::structures::contour_circulator;

typedef std::vector<point_type>::iterator point_iter;

namespace geom {
namespace predicates {

    using geom::structures::point_type;

    enum turn_type
    {
        COLLINEAR = 0, LEFT, RIGHT
    };

    namespace
    {
        template<typename Scalar>
        int sgn(Scalar x)
        {
            if (x == 0)
                return 0;
            else
                return x < 0 ? -1 : 1;
        }
    }

    turn_type turn( point_type const & a,
                    point_type const & b,
                    point_type const & c )
    {
        auto v1 = b - a;
        auto v2 = c - a;

        auto sign = sgn(v1 ^ v2);

        switch (sign)
        {
            case -1: return RIGHT;
            case  1: return LEFT;
            default: return COLLINEAR;
        }
    }
}}

namespace geom {
namespace structures {

    contour_type::contour_type(std::vector<point_type> && pts)
        : pts_(std::move(pts))
    {}

    struct contour_builder_type
    {
        void add_point(point_type const & pt)
        {
            pts_.push_back(pt);
        }

        contour_type get_result()
        {
            return contour_type(std::move(pts_));
        }

    private:
        std::vector<point_type> pts_;
    };
}}

namespace geom {
namespace algorithms {
namespace convex_hull {

contour_type merge_hull_impl(std::vector<point_type>::const_iterator beg,
                             std::vector<point_type>::const_iterator end);
contour_type merge(contour_type const & A, contour_type const & B);

using geom::structures::contour_builder_type;

contour_type merge_hull(std::vector<point_type> pts)
{
    if (pts.size() < 2) {
        throw std::logic_error("not enough points to build convex hull");
    }
    boost::sort(pts);

    return merge_hull_impl(pts.begin(), pts.end());
}

contour_type merge_hull_impl(std::vector<point_type>::const_iterator beg,
                             std::vector<point_type>::const_iterator end)
{
    size_t distance = std::distance(beg, end);
    if (distance <= 2) {
        contour_builder_type builder;
        for (std::vector<point_type>::const_iterator it = beg; it != end; ++it) {
            builder.add_point(*it);
        }
        return builder.get_result();
    }

    std::vector<point_type>::const_iterator middle = beg + distance / 2;
    contour_type left_contour = merge_hull_impl(beg, middle);
    contour_type right_contour = merge_hull_impl(middle, end);

    return merge(left_contour, right_contour);
}

bool is_tangent(point_type const & p1, point_type const & p2, contour_type const & contour)
{
    using namespace geom::predicates;

    for (auto it = contour.begin(); it != contour.end(); ++it) {
        if (turn(p1, p2, *it) == RIGHT) {
            return false;
        }
    }

    return true;
}

bool is_tangent(point_type const & p1, point_type const & p2, point_type const & p3)
{
    using namespace geom::predicates;

    return turn(p1, p2, p3) != RIGHT;
}

void find_tangent(contour_circulator & circular_A, contour_circulator & circular_B)
{
    while (!is_tangent(*circular_A, *circular_B, circular_A.prev())
           || !is_tangent(*circular_A, *circular_B, circular_B.next())) {
        while (!is_tangent(*circular_A, *circular_B, circular_A.prev())) {
            --circular_A;
        }
        while (!is_tangent(*circular_A, *circular_B, circular_B.next())) {
            ++circular_B;
        }
    }
}

void set_leftmost(contour_circulator & current)
{
    contour_circulator prev(current), next(current);
    --prev;
    ++next;

    if (prev == current && current == next) {
        return;
    }

    while (!(*current < *prev && *current < *next)) {
        ++current;
        ++prev;
        ++next;
    }
}

void set_rightmost(contour_circulator & current)
{
    contour_circulator prev(current), next(current);
    --prev;
    ++next;

    if (prev == current && current == next) {
        return;
    }

    while (!(*current > *prev && *current > *next)) {
        ++current;
        ++prev;
        ++next;
    }
}

bool is_on_one_line(contour_type const & contour) {
    using namespace geom::predicates;

    if (contour.vertices_num() <= 2) {
        return true;
    }

    for (auto it = contour.begin(); it != contour.end() - 2; ++it) {
        if (turn(*it, *(it + 1), *(it + 2)) != COLLINEAR) {
            return false;
        }
    }

    return true;
}

bool is_on_one_line(contour_type const & A, contour_type const & B)
{
    using namespace geom::predicates;

    contour_builder_type builder;
    for (auto it = A.begin(); it != A.end(); ++it) {
        builder.add_point(*it);
    }
    for (auto it = B.begin(); it != B.end(); ++it) {
        builder.add_point(*it);
    }
    contour_type contour = builder.get_result();

    return is_on_one_line(contour);
}

contour_type build_oneline_contour(contour_type const & A, contour_type const & B)
{
    contour_circulator circular_A(A);
    contour_circulator circular_B(B);
    set_leftmost(circular_A);
    set_leftmost(circular_B);

    contour_circulator end_A(circular_A);
    contour_circulator end_B(circular_B);
    --end_A;
    --end_B;

    contour_builder_type builder;
    while (circular_A != end_A) {
        builder.add_point(*circular_A);
        ++circular_A;
    }
    builder.add_point(*end_A);

    while (circular_B != end_B) {
        builder.add_point(*circular_B);
        ++circular_B;
    }
    builder.add_point(*end_B);

    return builder.get_result();
}

void add_oneline_contour_to_builder(contour_builder_type & builder,
                                    contour_circulator beg, contour_circulator end)
{
    if (++beg == end) {
        --beg;
        while (beg != end) {
            builder.add_point(*beg);
            --beg;
        }
    }
    else {
        --beg;
        while (beg != end) {
            builder.add_point(*beg);
            ++beg;
        }
    }
    builder.add_point(*end);
}

void add_convex_contour_to_builder(contour_builder_type & builder,
                                   contour_circulator beg, contour_circulator end)
{
    while (beg != end) {
        builder.add_point(*beg);
        ++beg;
    }
    builder.add_point(*end);
}

void add_contour_to_builder(contour_builder_type & builder, contour_type const & contour,
                            contour_circulator const & beg, contour_circulator const & end)
{
    if (is_on_one_line(contour)) {
        add_oneline_contour_to_builder(builder, beg, end);
    }
    else {
        add_convex_contour_to_builder(builder, beg, end);
    }
}

contour_type merge(contour_type const & A, contour_type const & B)
{
    if (is_on_one_line(A, B)) {
        return build_oneline_contour(A, B);
    }

    contour_circulator circular_A(A);
    contour_circulator circular_B(B);

    /*find lower tangent*/
    set_rightmost(circular_A);
    set_leftmost(circular_B);
    find_tangent(circular_A, circular_B);
    contour_circulator a_down(circular_A);
    contour_circulator b_down(circular_B);

    /*find upper tangent*/
    set_rightmost(circular_A);
    set_leftmost(circular_B);
    find_tangent(circular_B, circular_A);
    contour_circulator a_up(circular_A);
    contour_circulator b_up(circular_B);

    contour_builder_type builder;
    add_contour_to_builder(builder, A, a_up, a_down);
    add_contour_to_builder(builder, B, b_down, b_up);

    return builder.get_result();
}

}}}
