/*
 * Neuromapp - exception.cpp, Copyright (c), 2015,
 * Timothee Ewart - Swiss Federal Institute of technology in Lausanne,
 * timothee.ewart@epfl.ch,
 * All rights reserved.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library.
 */

/**
 * @file neuromapp/test/app/exception.cpp
 *  Test the generation of the exeption for C/C++ interface
 */

#define BOOST_TEST_MODULE ExceptionTest
#include <boost/test/unit_test.hpp>

#include "utils/error.h"
#include "app/driver_exception.h"

void helper_c_exception(int error){
    if(error)
        throw mapp::driver_exception(error);
}

void helper_cpp_exception(std::string const &  m){
        throw mapp::driver_exception(m);
}

BOOST_AUTO_TEST_CASE(exception){
    // exception generated by C interface
    BOOST_CHECK_THROW(helper_c_exception(mapp::MAPP_BAD_ARG), mapp::driver_exception);
    BOOST_CHECK_THROW(helper_c_exception(mapp::MAPP_USAGE), mapp::driver_exception);
    BOOST_CHECK_THROW(helper_c_exception(mapp::MAPP_BAD_DATA), mapp::driver_exception);
    BOOST_CHECK_THROW(helper_c_exception(mapp::MAPP_BAD_THREAD), mapp::driver_exception);
    BOOST_CHECK_THROW(helper_c_exception(666), mapp::driver_exception);
    // exception generated by C++ interface
    std::string message("Does it work?");
    BOOST_CHECK_THROW(helper_cpp_exception(message), mapp::driver_exception);
    try{
        helper_cpp_exception(message);
    }catch(mapp::driver_exception & e){
        BOOST_CHECK(e.what() == message);
    }
}