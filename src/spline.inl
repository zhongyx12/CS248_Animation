// Given a time between 0 and 1, evaluates a cubic polynomial with
// the given endpoint and tangent values at the beginning (0) and
// end (1) of the interval.  Optionally, one can request a derivative
// of the spline (0=no derivative, 1=first derivative, 2=2nd derivative).
template <class T>
inline T Spline<T>::cubicSplineUnitInterval(
    const T& position0, const T& position1, const T& tangent0,
    const T& tangent1, double normalizedTime, int derivative) {
  // TODO (Animation) Task 1a
  T result;
  double t3, t2, t1, t0;
  T a, b, c, d;

  t0 = 1.0;
  t1 = normalizedTime;
  t2 = t1 * t1;
  t3 = t1 * t2;

  if (derivative == 0) {
    a = 2 * position0 + tangent0 - 2 * position1 + tangent1;
    b = (-3) * position0 - 2 * tangent0 + 3 * position1 - 1 * tangent1;
    c = tangent0;
    d = position0;
    return (a * t3 + b * t2 + c * t1 + d * t0);
  }
  else if (derivative == 1) {
    b = 6 * position0 + 3 * tangent0 - 6 * position1 + 3 * tangent1;
    c = (-6) * position0 - 4 * tangent0 + 6 * position1 - 2 * tangent1;
    d = tangent0;
    return (b * t2 + c * t1 + d * t0);
  }
  else if (derivative == 2) {
    c = 12 * position0 + 6 * tangent0 - 12 * position1 + 6 * tangent1;
    d = (-6) * position0 - 4 * tangent0 + 6 * position1 - 2 * tangent1;
    return (c * t1 + d * t0);
  }
  else {
    return T();
  }
}

// Returns a state interpolated between the values directly before and after the
// given time.
template <class T>
inline T Spline<T>::evaluate(double time, int derivative) {
  // TODO (Animation) Task 1b
  if (knots.size() < 1) {
    return T();
  }
  else if (knots.size() == 1) {
    if (derivative == 0) {
      return knots.begin()->second;
    }
    else if (derivative == 1 || derivative == 2) {
      return 0 * T();
    }
    else {
      return T();
    }
  }
  else if (time <= knots.begin()->first) {
    if (derivative == 0) {
      return knots.begin()->second;
    }
    else if (derivative == 1 || derivative == 2) {
      return 0 * T();
    }
    else {
      return T();
    }
  }
  else if (time >= knots.rbegin()->first) {
    if (derivative == 0) {
      return knots.rbegin()->second;
    }
    else if (derivative == 1 || derivative == 2) {
      return 0 * T();
    }
    else {
      return T();
    }
  }

  auto iter2 = knots.upper_bound(time);
  auto iter1 = std::prev(iter2);
  double t1 = iter1->first, t2 = iter2->first;
  T p1 = iter1->second, p2 = iter2->second;
  double t0, t3;
  T p0, p3, m1, m2;
  if (iter1 == knots.begin()) {
    t0 = t1 - (t2 - t1);
    p0 = p1 - (p2 - p1);
  } else {
    auto iter0 = std::prev(iter1);
    t0 = iter0->first;
    p0 = iter0->second;
  }
  if (std::next(iter2) == knots.end()) {
    t3 = t2 + (t2 - t1);
    p3 = p2 + (p2 - p1);
  } else {
    auto iter3 = std::next(iter2);
    t3 = iter3->first;
    p3 = iter3->second;
  }
  m1 = (p2 - p0) / (t2 - t0);
  m2 = (p3 - p1) / (t3 - t1);

  double L = t2 - t1;
  double _time = (time - t1) / L;
  T _m1 = m1 * L;
  T _m2 = m2 * L;
  T tmp_result = Spline<T>::cubicSplineUnitInterval(p1, p2, _m1, _m2, _time, derivative);

  if (derivative == 0) {
    return tmp_result;
  }
  else if (derivative == 1) {
    return tmp_result / L;
  }
  else if (derivative == 2) {
    return tmp_result / (L * L);
  }
  else {
    return tmp_result;
  }
}

// Removes the knot closest to the given time,
//    within the given tolerance..
// returns true iff a knot was removed.
template <class T>
inline bool Spline<T>::removeKnot(double time, double tolerance) {
  // Empty maps have no knots.
  if (knots.size() < 1) {
    return false;
  }

  // Look up the first element > or = to time.
  typename std::map<double, T>::iterator t2_iter = knots.lower_bound(time);
  typename std::map<double, T>::iterator t1_iter;
  t1_iter = t2_iter;
  t1_iter--;

  if (t2_iter == knots.end()) {
    t2_iter = t1_iter;
  }

  // Handle tolerance bounds,
  // because we are working with floating point numbers.
  double t1 = (*t1_iter).first;
  double t2 = (*t2_iter).first;

  double d1 = fabs(t1 - time);
  double d2 = fabs(t2 - time);

  if (d1 < tolerance && d1 < d2) {
    knots.erase(t1_iter);
    return true;
  }

  if (d2 < tolerance && d2 < d1) {
    knots.erase(t2_iter);
    return t2;
  }

  return false;
}

// Sets the value of the spline at a given time (i.e., knot),
// creating a new knot at this time if necessary.
template <class T>
inline void Spline<T>::setValue(double time, T value) {
  knots[time] = value;
}

template <class T>
inline T Spline<T>::operator()(double time) {
  return evaluate(time);
}
