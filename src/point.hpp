// ============================================================================
// Copyright (C) SlideforMAP++. All rights reserved.
//
// Authors: Denis Cohen-Corticchiato (DOHCC)
// Email:   denis.cohen@gmail.com
//
// This file is part of SlideforMAP++  
// ============================================================================

#ifndef POINT_HPP
#define POINT_HPP

#include <iostream>
#include <cmath>

//=============================================================================

class Point
{
   public:

      Point() { };
      Point(double x, double y): x_(x), y_(y) { };
      Point(double x, double y, double z): x_(x), y_(y), z_(z) { };

      //
      // Operations
      //
      double dist2(const Point & v) const
      {
         auto dx = this->x_ - v.x_;
         auto dy = this->y_ - v.y_;
         return dx * dx + dy * dy;
      }

      double dist(const Point & v) const
      {
         return sqrt(dist2(v));
      }

      double dist(double x, double y) const
      {
         double dx = this->x_ - x;
         double dy = this->y_ - y;
         return dx * dx + dy * dy;
      }

      double norm2() const
      {
         return x_ * x_ + y_ * y_;
      }

      double length() const
      {
         return sqrt(x_ * x_ + y_ * y_);
      }

      void set(double x, double y)
      {
         this->x_ = x;
         this->y_ = y;
      }

      void add(double x, double y)
      {
         this->x_ += x;
         this->y_ += y;
      }

      double x_{0.0};
      double y_{0.0};
      double z_{0.0};
};

inline bool operator < (const Point& v1, const Point& v2)
{
   return v2.length() > v1.length();
}

inline bool operator == (const Point& v1, const Point& v2)
{
   return (v1.x_ == v2.x_) && (v1.y_ == v2.y_);
}

//inline bool almost_equal(const Point& v1, const Point& v2, int ulp=2)
//{
//   return almost_equal(v1.x_, v2.x_, ulp) && almost_equal(v1.y_, v2.y_, ulp);
//}

#endif
