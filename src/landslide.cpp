// ============================================================================
// Copyright (C) SlideforMAP++. All rights reserved.
//
// Authors: Denis Cohen-Corticchiato (DOHCC)
// Email:   denis.cohen@gmail.com
//
// This file is part of SlideforMAP++  
// ============================================================================

#include "landslide.hpp"

Landslide::Landslide()
{
}

Landslide::~Landslide()
{
}

double Landslide::GetSlopeDeg() {return slope_ * 180 / M_PI;}

void Landslide::SetFrictionAngle(const unsigned int & angle)
{
	friction_angle_ = M_PI*angle/180;
    while (friction_angle_ < 0) friction_angle_+= 2*M_PI;
    while (friction_angle_ > M_PI) friction_angle_ -= M_PI;
    if (friction_angle_ == M_PI/2) friction_angle_ = 179*M_PI/180;
    if (friction_angle_ > M_PI/2) friction_angle_ = M_PI - friction_angle_;
    if (friction_angle_ == 0) friction_angle_ = M_PI/180;
}

