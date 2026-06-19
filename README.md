# star_ship_fly

## Project Structure

```
project/
├── src/
│   ├── imu_fx_10hz.c
│   ├── imu_gy_10hz.c
│   ├── kalman_filter.c
│   ├── star_ship_main.c
│   ├── ins_01.c
│   ├── ins_02.c
│   ├── matrix_proc.c
│   └── navigation_process_function.c
├── inc/
│   └── *.h
├── build.sh
└── README.md
```

## Build Instructions

### Linux / macOS / Git Bash (Windows)

```bash
chmod +x build.sh
./build.sh
```

### Windows (Command Prompt)

```bat
build.bat
```

## Compiler Flags

| Flag | Description |
|------|-------------|
| `src/*.c` | Source files located in the `src/` folder |
| `-I inc` | Header files located in the `inc/` folder |
| `-o star_ship_fly.exe` | Output executable name |

## Requirements

- GCC compiler
- Windows users: install [Git for Windows](https://git-scm.com/download/win) to run `.sh` scripts