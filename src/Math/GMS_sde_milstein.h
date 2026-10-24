


/*
Copyright © 2011-2017, Andrew D. Horchler
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of Case Western Reserve University nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL ANDREW D. HORCHLER BE LIABLE FOR ANY
   DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/* Copyright (C) Bernard Gingold, 2020-2026 
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*
   Ported to the C++17 by Bernard Gingold.
*/
#include <Eigen/Dense>
#include <functional>
#include <vector>
#include <stdexcept>
#include <random>
#include <cmath>

namespace gms
{

namespace math 
{
    struct SdeOptions
    {
        // Matches the original solver defaults:
        // default is Stratonovich, not Ito
        bool stratonovich = true;

        // Original solver supports diagonal and correlated noise.
        // This faithful C++17 core implementation supports diagonal noise.
        bool diagonal_noise = true;

        // Constant drift/diffusion detection
        bool const_ffun = false;
        bool const_gfun = false;
        bool const_dgfun = false;

        // Optional callback functions
        std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)> f;
        std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)> g;
        std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)> dg;
    };

    struct SdeMilsteinResult
    {
        Eigen::MatrixXd Y;  // solution, size lt x N
        Eigen::MatrixXd W;  // integrated Wiener increments, size lt x D
    };

    inline Eigen::VectorXd eval_vector_or_constant(
        const std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)>& fn,
        double t,
        const Eigen::VectorXd& y,
        const Eigen::VectorXd& constant,
        bool is_constant)
    {
        if (is_constant)
            return constant;
        return fn(t, y);
    }

    // Faithful central solver from sde_milstein.m
    inline SdeMilsteinResult sde_milstein(
        const std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)>& f,
        const std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)>& g,
        const std::vector<double>& tspan,
        const Eigen::VectorXd& y0,
        const SdeOptions& options = {})
    {
        if (tspan.size() < 2)
            throw std::invalid_argument("tspan must contain at least two time points");

        const std::size_t lt = tspan.size();
        const std::size_t N = static_cast<std::size_t>(y0.size());

        // In the original file, the general solver supports N-dimensional diagonal noise.
        // We mirror the diagonal-noise case directly.
        const std::size_t D = N;

        Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(static_cast<int>(lt), static_cast<int>(N));
        Eigen::MatrixXd W = Eigen::MatrixXd::Zero(static_cast<int>(lt), static_cast<int>(D));

        Y.row(0) = y0.transpose();

        std::mt19937 rng( std::random_device{}() );
        std::normal_distribution<double> normal(0.0, 1.0);

        double Ti = tspan.front();
        Eigen::VectorXd Yi = y0;

        // We keep the original numerical structure of the solver:
        // - compute dt
        // - draw Wiener increments dW
        // - evaluate f(Ti,Yi)
        // - evaluate g(Ti,Yi)
        // - use derivative-free or derivative-based Milstein update
        for (std::size_t i = 0; i + 1 < lt; ++i)
        {
            const double dt = tspan[i + 1] - tspan[i];
            const double sdt = std::sqrt(std::abs(dt));

            // Wiener increment vector dW (diagonal noise)
            Eigen::VectorXd dW = Eigen::VectorXd::Zero(static_cast<int>(D));
            for (std::size_t j = 0; j < D; ++j)
                dW[j] = sdt * normal(rng);

            // integrated Wiener increments: W(i+1,:) = W(i,:) + dW
            W.row(static_cast<int>(i + 1)) = W.row(static_cast<int>(i)) + dW.transpose();

            // Evaluate drift and diffusion
            Eigen::VectorXd fout = f(Ti, Yi) * dt;
            Eigen::VectorXd gout = g(Ti, Yi);

            Eigen::VectorXd dW2;
            Eigen::VectorXd Ybar;

            if (options.dg)
            {
                // Explicit Milstein with derivative of g
                // Original logic:
                // if Stratonovich
                //     dW2 = 0.5 * dW.^2;
                // else
                //     dW2 = 0.5 * (dW.^2 - dt);
                //
                // Yi = Yi + fout + gout .* (dW + dgout .* dW2);

                Eigen::VectorXd dgout = options.dg(Ti, Yi);

                if (options.stratonovich)
                    dW2 = 0.5 * dW.array().square();
                else
                    dW2 = 0.5 * (dW.array().square() - dt);

                Yi = Yi + fout + gout.array() * (dW.array() + dgout.array() * dW2.array());
            }
            else
            {
                // Derivative-free Milstein
                // Original logic:
                // if Stratonovich
                //     dW2 = (0.5/sdt) * dW.^2;
                // else
                //     dW2 = (0.5/sdt) * (dW.^2 - dt);
                //
                // Ybar = Yi + gout .* dW;
                // Yi = Ybar + fout + (g(Ti, Ybar) - gout) .* dW2;

                if (options.stratonovich)
                    dW2 = (0.5 / sdt) * dW.array().square();
                else
                    dW2 = (0.5 / sdt) * (dW.array().square() - dt);

                Ybar = Yi + gout.array() * dW.array();

                Eigen::VectorXd gbar = g(Ti, Ybar);
                Yi = Ybar + fout + (gbar - gout).array() * dW2.array();
            }

            // Store current solution
            Y.row(static_cast<int>(i + 1)) = Yi.transpose();

            // Advance time
            Ti = tspan[i + 1];
        }

        return {Y, W};
    }

    // Convenience wrapper to mimic the original function structure:
    // SDE_MILSTEIN(f, g, tspan, y0, options)
    inline SdeMilsteinResult sde_milstein(
        const std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)>& f,
        const std::function<Eigen::VectorXd(double, const Eigen::VectorXd&)>& g,
        const std::vector<double>& tspan,
        const Eigen::VectorXd& y0,
        bool stratonovich)
    {
        SdeOptions options;
        options.stratonovich = stratonovich;
        options.f = f;
        options.g = g;
        return sde_milstein(f, g, tspan, y0, options);
    }
    

#include <iostream>


int unit_test()
{
    using namespace sdetools;

    const double mu = 1.0;
    Eigen::Vector2d sig;
    sig << 0.1, 0.5;

    std::vector<double> t;
    const double dt = 1e-2;
    for (double tt = 0.0; tt <= 1.0 + 1e-12; tt += dt)
        t.push_back(tt);

    auto f = [mu](double, const Eigen::VectorXd& y) -> Eigen::VectorXd
    {
        return mu * y;
    };

    auto g = [sig](double, const Eigen::VectorXd& y) -> Eigen::VectorXd
    {
        return sig.array() * y.array();
    };

    // Derivative of noise term: dg = diag(sig)
    auto dg = [sig](double, const Eigen::VectorXd&) -> Eigen::VectorXd
    {
        // This is a diagonal noise derivative in the vectorized sense:
        // d/dy [g(t,y)] = diag(sig)
        return sig;
    };

    // Original MATLAB example uses DGFUN to trigger the explicit Milstein method
    SdeOptions options;
    options.stratonovich = true;
    options.diagonal_noise = true;
    options.dg = dg;

    auto result = sde_milstein(f, g, t, Eigen::Vector2d(1.0, 1.0), options);

    std::cout << "Y:\n" << result.Y << "\n\n";
    std::cout << "W:\n" << result.W << "\n";

    return 0;
}

}  // namespace math
} // namespace gms
