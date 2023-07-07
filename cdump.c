#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

struct atom {
    uint64_t id;
    uint64_t type;
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
};

void print_atom(const struct atom * const a)
{
    printf("%" PRIu64 ", %" PRIu64 ", %.17g, %.17g, %.17g, %.17g, %.17g, %.17g\n",
            a->id, a->type, a->x, a->y, a->z, a->vx, a->vy, a->vz);
}

enum parser_state {
    WAIT_FOR_ITEM = 0,
    IN_TIMESTEP = 1,
    IN_NUM_ATOMS = 2,
    IN_BOX_BOUNDS = 3,
    IN_ATOMS = 4,
};

struct dump_parser {
    uint64_t timestep;
    uint64_t num_atoms;
    uint64_t num_bounds;
    uint64_t read_atoms;
    enum parser_state state;
};

const struct dump_parser DEFAULT_DUMP_PARSER = {
    .timestep = 0,
    .num_atoms = 0,
    .num_bounds = 0,
    .read_atoms = 0,
    .state = WAIT_FOR_ITEM,
};

struct fixed_string {
    const char *s;
    size_t n;
};

#define MAKE_FIXED_STRING(s) (struct fixed_string){s, sizeof(s) - 1}
#define DIM(x) (sizeof(x) / sizeof(x[0]))

enum parser_state next_state_from_known_line(const char *line, size_t line_length)
{
    // printf("LINE: %d %.*s\n", (int)line_length, (int)line_length, line);
    // printf("RLINE: %d %s\n", (int)sizeof("ITEM: TIMESTEP"), "ITEM: TIMESTEP");

    const struct fixed_string known[] = {
        [WAIT_FOR_ITEM] = MAKE_FIXED_STRING(""),
        [IN_TIMESTEP]   = MAKE_FIXED_STRING("ITEM: TIMESTEP"),
        [IN_NUM_ATOMS]  = MAKE_FIXED_STRING("ITEM: NUMBER OF ATOMS"),
        [IN_BOX_BOUNDS] = MAKE_FIXED_STRING("ITEM: BOX BOUNDS pp pp ss"),
        [IN_ATOMS]      = MAKE_FIXED_STRING("ITEM: ATOMS id type x y z vx vy vz"),
    };

    // prefix search to ignore junk at end.
    for (size_t i = 1; i < DIM(known); ++i) {
        if (line_length >= known[i].n) {
            if (memcmp(line, known[i].s, known[i].n) == 0) {
                return (enum parser_state)i;
            }
        }
    }
    return WAIT_FOR_ITEM;
}

bool parse_atom(struct atom *a, const char *line, size_t line_length)
{
    assert(a != NULL);
    assert(line != NULL);

    int last_space = 0;
    int which_token = 0;
    for (size_t i = 0; i < line_length; ++i) {
        if (line[i] == ' ' || line[i] == '\n' || line[i] == '\0') {
            switch (which_token) {
                case 0: {
                    a->id = strtoull(&line[last_space], NULL, 10);
                } break;
                case 1: {
                    a->type = strtoull(&line[last_space], NULL, 10);
                } break;
                case 2: {
                    a->x= strtod(&line[last_space], NULL);
                } break;
                case 3: {
                    a->y= strtod(&line[last_space], NULL);
                } break;
                case 4: {
                    a->z= strtod(&line[last_space], NULL);
                } break;
                case 5: {
                    a->vx= strtod(&line[last_space], NULL);
                } break;
                case 6: {
                    a->vy= strtod(&line[last_space], NULL);
                } break;
                case 7: {
                    a->vz= strtod(&line[last_space], NULL);
                } break;
            }

            last_space = i + 1;
            which_token++;

            if (which_token >= 8) {
                return true;
            }
        }
    }
    return false;
}

static struct atom atoms[1048576] = {0};

double cross_sectional_area(double radius, double zsphere, double zslice)
{
    const double dh = fabs(zsphere - zslice);
    if (dh >= radius) {
        return 0.0;
    }

    const double a = sqrt(radius * radius - dh * dh);
    return M_PI * a * a;
}

double dvx2_over_slice(struct atom *atoms, size_t num_atoms, double grain_radius, double zslice)
{
    double dvx2 = 0;
    double A_slice = 0;
    double vx_slice = 0;
    for (size_t i = 0; i < num_atoms; ++i) {
        const struct atom * const entry = &atoms[i];
        const double A = cross_sectional_area(grain_radius, entry->z, zslice);
        vx_slice += A * entry->vx;
        A_slice += A;
    }

    if (A_slice <= 0) {
        return 0;
    }
    vx_slice /= A_slice;

    for (size_t i = 0; i < num_atoms; ++i) {
        const struct atom * const entry = &atoms[i];
        const double A = cross_sectional_area(grain_radius, entry->z, zslice);
        const double dvx = (entry->vx - vx_slice);
        dvx2 += A * dvx * dvx;
    }

    dvx2 /= A_slice;
    return dvx2;
}

double dv2_over_slice(struct atom *atoms, size_t num_atoms, double grain_radius, double zslice)
{
    double dvx2 = 0;
    double dvy2 = 0;
    double dvz2 = 0;
    double A_slice = 0;
    double vx_slice = 0;
    double vy_slice = 0;
    double vz_slice = 0;
    for (size_t i = 0; i < num_atoms; ++i) {
        const struct atom * const entry = &atoms[i];
        const double A = cross_sectional_area(grain_radius, entry->z, zslice);
        vx_slice += A * entry->vx;
        vy_slice += A * entry->vy;
        vz_slice += A * entry->vz;
        A_slice += A;
    }

    if (A_slice <= 0) {
        return 0;
    }
    vx_slice /= A_slice;
    vy_slice /= A_slice;
    vz_slice /= A_slice;

    for (size_t i = 0; i < num_atoms; ++i) {
        const struct atom * const entry = &atoms[i];
        const double A = cross_sectional_area(grain_radius, entry->z, zslice);
        const double dvx = (entry->vx - vx_slice);
        const double dvy = (entry->vy - vy_slice);
        const double dvz = (entry->vz - vz_slice);
        dvx2 += A * dvx * dvx;
        dvy2 += A * dvy * dvy;
        dvz2 += A * dvz * dvz;
    }

    dvx2 /= A_slice;
    dvy2 /= A_slice;
    dvz2 /= A_slice;
    return (dvx2 + dvy2 + dvz2);
}


void process_frame(struct atom *atoms, size_t num_atoms, FILE *fout)
{

/*
    const double grain_radius = 0.5;
    double dvx2_at_slices[100] = {0};
    const double slice_dh = 30.0 * grain_radius / DIM(dvx2_at_slices);
    for (size_t i = 0; i < DIM(dvx2_at_slices); ++i) {
        const double zslice = i * slice_dh;
        dvx2_at_slices[i] = dvx2_over_slice(atoms, num_atoms, grain_radius, zslice);
    }

    for (size_t i = 0; i < DIM(dvx2_at_slices); ++i) {
        fprintf(fout, "%.17g,", dvx2_at_slices[i]);
    }
    fprintf(fout, "\n");
*/

    const double grain_radius = 0.5;
    double dv2_at_slices[100] = {0};
    const double slice_dh = 30.0 * grain_radius / DIM(dv2_at_slices);
    for (size_t i = 0; i < DIM(dv2_at_slices); ++i) {
        const double zslice = i * slice_dh;
        dv2_at_slices[i] = dv2_over_slice(atoms, num_atoms, grain_radius, zslice);
    }

    for (size_t i = 0; i < DIM(dv2_at_slices); ++i) {
        fprintf(fout, "%.17g,", dv2_at_slices[i]);
    }
    fprintf(fout, "\n");
}

int main(int argc, char *argv[])
{
    assert(argc >= 2);
    FILE *fin = fopen(argv[1], "r");
    assert(fin != NULL);

    FILE *fout = stdout;
    if (argc >= 3) {
        fout = fopen(argv[2], "w+");
        assert(fout != NULL);
    }

    struct dump_parser s = DEFAULT_DUMP_PARSER;
    size_t remaining = 0;
    char line[2048] = {0};
    size_t num_frames = 0;
    while (!feof(fin)) {
        assert(remaining < sizeof(line));
        size_t num_read = fread(line + remaining, 1, sizeof(line) - 1 - remaining, fin);
        // printf("%c", line[0]);
        // sum += line[0];

        // printf("READ BYTES: %zu\n", num_read);

        remaining += num_read;
        // printf("REMAINING: %zu\n", remaining);
        // printf("AFTER READ: |%.*s|\n", (int)remaining, line);
        // for (size_t i = 0; i < remaining; ++i) {
        //     printf("%02X ", line[i]);
        // }
        // printf("\n");
        while (remaining > 0) {
            size_t line_end = sizeof(line);
            for (size_t i = 0; i < remaining; ++i) {
                if (line[i] == '\n' || line[i] == '\0') {
                    line_end = i;
                    break;
                }
            }

            if (line_end == sizeof(line)) {
                // printf("NO END IN LINE\n");
                break;
            }

            // printf("LINE TO PROCESS: |%.*s|\n", (int)line_end, line);

            switch (s.state) {
                case WAIT_FOR_ITEM: {
                    s.state = next_state_from_known_line(line, line_end);
                } break;
                case IN_TIMESTEP: {
                    char z[64] = {0};
                    assert(line_end < sizeof(z));
                    memcpy(z, line, line_end);
                    s.timestep = strtoull(z, NULL, 10);
                    s.state = WAIT_FOR_ITEM;
                } break;
                case IN_NUM_ATOMS: {
                    char z[64] = {0};
                    assert(line_end < sizeof(z));
                    memcpy(z, line, line_end);
                    s.num_atoms = strtoull(z, NULL, 10);
                    s.state = WAIT_FOR_ITEM;
                    // printf("num atoms: %zu\n", s.num_atoms);
                } break;
                case IN_BOX_BOUNDS: {
                    s.num_bounds++;
                    if (s.num_bounds >= 3) {
                        s.state = WAIT_FOR_ITEM;
                    }
                } break;
                case IN_ATOMS: {
                    // printf("read atoms: %zu\n", s.read_atoms);
                    assert(s.read_atoms <= DIM(atoms));
                    if (parse_atom(&atoms[s.read_atoms], line, line_end)) {
                        // print_atom(&atoms[s.read_atoms]);
                        s.read_atoms++;
                    } else {
                        printf("atom failed to parse: %.*s", (int)line_end, line);
                    }

                    if (s.num_atoms == s.read_atoms) {
                        process_frame(atoms, s.num_atoms, fout);
                        // printf("process frame %zu\n", num_frames);
                        num_frames++;
                        s = DEFAULT_DUMP_PARSER;
                    }
                } break;
                default:
                    printf("unknown state");
                    assert(false);
            }

            assert(remaining > line_end);
            remaining -= line_end;
            assert(remaining >= 1);
            remaining--;
            assert(line_end < sizeof(line));
            // printf("MOVE BYTES: %zu\n", remaining);
            // printf("MOVE TO FRONT: |%.*s|\n", (int)remaining, &line[line_end + 1]);

            memmove(line, &line[line_end + 1], remaining);

            // printf("REMAINING: %zu\n", remaining);
        }
    }

    if (fout != stdout) {
        fclose(fout);
    }

    return 0;
}
