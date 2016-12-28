//===========================================================================
//
// File: StopWatch.cpp
//
// Created: Thu Jul  2 23:04:51 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#if HAVE_CONFIG_H
#include "config.h"
#endif
#include <boost/date_time/posix_time/posix_time.hpp>
#include <opm/core/utility/StopWatch.hpp>
#include <opm/common/ErrorMacros.hpp>

namespace Opm
{

    namespace time
    {
        using Duration = boost::posix_time::time_duration;

	StopWatch::StopWatch()
	    : state_(State::UnStarted)
	{
	}


	void StopWatch::start()
	{
	    start_time_ = currentTime();
	    last_time_ = start_time_;
	    state_ = State::Running;
	}

	void StopWatch::stop()
	{
	    if (state_ != State::Running) {
		OPM_THROW(std::runtime_error, "Called stop() on a StopWatch that was not running.");
	    }
	    stop_time_ = currentTime();
	    state_ = State::Stopped;
	}

	double StopWatch::secsSinceLast()
	{
	    TimePoint run_time;
	    if (state_ == State::Running) {
		run_time = currentTime();
	    } else if (state_ == State::Stopped) {
		run_time = stop_time_;
	    } else {
		assert(state_ == State::UnStarted);
		OPM_THROW(std::runtime_error, "Called secsSinceLast() on a StopWatch that had not been started.");
	    }
	    Duration dur = run_time - last_time_;
	    last_time_ = run_time;
	    return double(dur.total_microseconds())/1000000.0;
	}

	double StopWatch::secsSinceStart()
	{
	    TimePoint run_time;
	    if (state_ == State::Running) {
		run_time = currentTime();
	    } else if (state_ == State::Stopped) {
		run_time = stop_time_;
	    } else {
		assert(state_ == State::UnStarted);
		OPM_THROW(std::runtime_error, "Called secsSinceStart() on a StopWatch that had not been started.");
	    }
	    Duration dur = run_time - start_time_;
	    return double(dur.total_microseconds())/1000000.0;
	}

        StopWatch::TimePoint StopWatch::currentTime() const
        {
            return boost::posix_time::microsec_clock::local_time();
        }

    } // namespace time

} // namespace Opm



