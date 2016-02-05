/*
 * Neuromapp - keyvalue.h, Copyright (c), 2015,
 * Judit Planas - Swiss Federal Institute of technology in Lausanne,
 * judit.planas@epfl.ch,
 * All rights reserved.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details..  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library.
 */

/**
 * @file neuromapp/keyvalue/keyvalue.h
 * keyvalue Miniapp
 */

#ifndef MAPP_KEYVALUE_EXECUTE_
#define MAPP_KEYVALUE_EXECUTE_

/** \fn keyvalue_execute(int argc, char *const argv[])
    \brief keyvalue Miniapp
    \param argc number of argument from the command line
    \param argv the command line from the driver or external call
    \return error message from mapp::mapp_error
 */
int keyvalue_execute(int argc, char* const argv[]);

#endif
