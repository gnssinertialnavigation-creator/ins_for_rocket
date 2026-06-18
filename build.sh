#!/bin/bash
 
gcc \
    src/rocket_body_force_10hz.c \
    src/rocket_body_omega_10hz.c \
    src/kalman_filter.c \
    src/kalman_filter_02.c \
    src/star_ship_main.c \
    src/ins_01.c \
    src/ins_02.c \
    src/ins_03.c \
    src/ins_04.c \
    src/matrix_proc.c \
    src/navigation_process_function.c \
    -I inc \
    -o star_ship_fly.exe