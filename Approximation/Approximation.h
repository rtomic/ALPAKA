/*
 *  Approximation.h
 *
 *  Created on: 24.02.2025
 *  Last changed: 11.06.2025
 *      Author: kempf
 *
 *      Class Approximation: Wrapper class for different approximation types
 *
 *
 *
 */

#pragma once

#include <vector>
#include <functional>



class Approximation
{
	public:
		virtual double at (std::vector<double> const &x) const = 0; // Wert der Approximation an Punkt x
		virtual double atUpToLevel(std::vector<double> const &x, std::size_t const &level) const = 0;
		virtual void solve(std::function<double(const std::vector<double> &)> const &targetFunc) = 0;
		virtual void solve(std::vector<double> const &targetFunc) = 0;
		virtual ~Approximation () = default; // Virtueller Destruktor für Polymorphismus

};
