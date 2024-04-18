
#include "../geometry/spline.h"

template<typename T> T Spline<T>::at(float time) const {

	// A4T1b: Evaluate a Catumull-Rom spline

	// Given a time, find the nearest positions & tangent values
	// defined by the control point map.

	// Transform them for use with cubic_unit_spline

	// Be wary of edge cases! What if time is before the first knot,
	// before the second knot, etc...

    // Step 1, knots <= 1
    // 1. no point
    if (!any()){
        return T();
    }
    // 2. only 1 point
    if (knots.size() == 1){
        return knots.begin()->second;
    }
    // 3. time less or equal than initial knot
    float initial_t = knots.begin()->first;
    if (time <= initial_t){
        return knots.begin()->second;
    }

    // 4. time greater or equal than initial knot
    auto last_knot_it = std::prev(knots.end());
    float end_t = last_knot_it->first;
    if (time >= end_t){
        return last_knot_it->second;
    }

    // Step2 knots >= 2
    // 0. time is the control point
    if (has(time)){
        return knots.at(time);
    }

    std::set<float> times = keys();
    std::vector<float> left_times, right_times;
    for (const float& t : times) {
        if (t < time) {
            left_times.push_back(t);
        } else if (t > time) {
            right_times.push_back(t);
        }
    }
    // 1. knots <4, mirror (which means either left or right or both will have no more than 2 points
    // 2. knots >= 4, call cubic_unit_spline
    // note: 2 can include 1 after 1 create virtual point

    float t0, t1, t2 , t3;
    T p0, p1, p2, p3;
    T m0, m1;

    t1 = left_times.back();
    t2 = right_times[0];
    p1 = knots.at(t1);
    p2 = knots.at(t2);
    // a) only 1 point to the left (instead of should be 2)
    if (left_times.size() <= 1){
        // create a virtual point
        t0 = t1 - (t2-t1);
        p0 = p1 - (p2-p1);
    }
    else {
        t0 = left_times[left_times.size() - 2];
        p0 = knots.at(t0);
    }
    // b) only 1 point to the right (instead of should be 2)
    if (right_times.size() <= 1){
        // create a virtual point
        t3 = t2 + (t2-t1);
        p3 = p2 + (p2-p1);
    }
    else {
        t3 = right_times[1];
        p3 = knots.at(t3);
    }
    m0 = (p2-p0)/(t2-t0);
    m1 = (p3-p1)/(t3-t1);
    // start point:p1 (not p0!), end point: p2 (not p3!)
    float t_prime = (time-t1)/(t2-t1);
    T m0_prime = m0 * (t2 - t1);
    T m1_prime = m1 * (t3 - t2);
    return cubic_unit_spline(t_prime, p1, p2, m0_prime, m1_prime);

//	return cubic_unit_spline(0.0f, T(), T(), T(), T());
}

template<typename T>
T Spline<T>::cubic_unit_spline(float time, const T& position0, const T& position1,
                               const T& tangent0, const T& tangent1) {

	// A4T1a: Hermite Curve over the unit interval

	// Given time in [0,1] compute the cubic spline coefficients and use them to compute
	// the interpolated value at time 'time' based on the positions & tangents

	// Note that Spline is parameterized on type T, which allows us to create splines over
	// any type that supports the * and + operators.
    float t = time;
    T p0 = position0;
    T p1 = position1;
    T m0 = tangent0;
    T m1 = tangent1;

    float h00 = 2*(float)pow(t, 3) - 3*(float)pow(t, 2) + 1;
    float h10 = (float)pow(t, 3) - 2*(float)pow(t, 2) + t;
    float h01 = -2*(float)pow(t, 3) + 3*(float)pow(t, 2);
    float h11 = (float)pow(t, 3) - (float)pow(t, 2);

    T pt = h00*p0 + h10*m0 + h01*p1 + h11*m1;

    return pt;
//	return T();
}

template class Spline<float>;
template class Spline<double>;
template class Spline<Vec4>;
template class Spline<Vec3>;
template class Spline<Vec2>;
template class Spline<Mat4>;
template class Spline<Spectrum>;
