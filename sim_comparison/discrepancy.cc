#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <set>
#include <vector>
#include <omp.h>

struct point {
  point(float e_, float i_, float u_, int idx_) : e(e_), i(i_), u(u_), idx(idx_) {}
  point() : e(0), i(0), u(0), idx(-1) {}
  point(const point& other) : e(other.e), i(other.i), u(other.u), idx(other.idx) {}
  float e, i, u;
  int idx;
};

bool get_col_int(const char* line_c, int col_no, int *val, const char *delim=" ") {
  char *line = strdup(line_c), *token;
  token = NULL;
  int col = 0;
  token = strtok(line, delim);
  while (token != NULL) {
    if (col++ == col_no) {
      *val = atoi(token);
      free(line);
      return true;
    }
    token = strtok(NULL, delim);
  }

  *val = atoi(token);
  free(line);
  return false;
}

bool get_col(const char* line_c, int col_no, float *val, const char *delim=" ") {
  char *line = strdup(line_c), *token;
  token = NULL;
  int col = 0;
  token = strtok(line, delim);
  while (token != NULL) {
    if (col++ == col_no) {
      *val = atof(token);
      free(line);
      return true;
    }
    token = strtok(NULL, delim);
  }

  *val = atof(token);
  free(line);
  return false;
}

std::vector<point> read_inf(const char* inf, int idxcol, int ecol, int icol, int ucol) {
  std::vector<point> ret;
  char *line = nullptr;
  size_t len = 0;
  ssize_t read;

  FILE* fp = fopen(inf, "r");
  if (fp == nullptr) {
    perror("Error in reading file");
    return ret;
  }

  int line_num = 0;
  while ((read = getline(&line, &len, fp)) != -1) {
    line_num++;
    int idx;
    if (!get_col_int(line, idxcol, &idx)) {
      fprintf(stderr, "Error on line %d: Could not find index column %d!! on line %s\n", line_num, idxcol, line);
      continue;
    }
    float en, ir, uir;
    if (!get_col(line, ecol, &en)) {
      fprintf(stderr, "Error on line %d: Could not find energy column %d!! on line %s\n", line_num, ecol, line);
      continue;
    }
    if (!get_col(line, icol, &ir)) {
      fprintf(stderr, "Error on line %d: Could not find irmsd column %d!! on line %s\n", line_num, icol, line);
      continue;
    }
    if (!get_col(line, ucol, &uir)) {
      fprintf(stderr, "Error on line %d: Could not find uirmsd column %d!! on line %s\n", line_num, ucol, line);
      continue;
    }
    // Ignore the input.
    if (idx == 9999) continue;

    ret.emplace_back(en, ir, uir, idx);
  }

  fclose(fp);
  if (line) free(line);

  return ret;
}

float vol(const point& bottom, const point& top, const point& sd) {
  return sqrt(pow((bottom.e - top.e) / sd.e, 2.0) +
      pow((bottom.i - top.i) / sd.i, 2.0) +
      pow((bottom.u - top.u) / sd.u, 2.0));
}

float one_disc(
    const std::vector<point> &pts,
    const point& zero,
    const point& corner, float max_vol, const point& sd) {
  int num_inside = 0;
  for (const auto& p : pts) {
    if (p.e <= corner.e && p.i <= corner.i && p.u <= corner.u) {
      num_inside++;
    }
  }

  // float first_part = (float)num_inside / pts.size(),
  //       second_part = vol(zero, corner, sd) / max_vol;
  // printf("d %f %f %f (%f,%f -> %f,%f)\n",
  //     first_part, second_part, std::abs(first_part - second_part),
  //     zero.e, zero.i, corner.e, corner.i);
  return std::abs((float)num_inside / pts.size() - vol(zero, corner, sd) / max_vol);
}

float dispersion(const std::vector<point> &pts) {
  // Get the min/max in each direction.
  point minp = pts[0], maxp = pts[0], meanp, stddevp;
  for (const auto& p : pts) {
    minp.e = std::min(p.e, minp.e);
    minp.i = std::min(p.i, minp.i);
    minp.u = std::min(p.u, minp.u);
    maxp.e = std::max(p.e, maxp.e);
    maxp.i = std::max(p.i, maxp.i);
    maxp.u = std::max(p.u, maxp.u);

    meanp.e += p.e;
    meanp.i += p.i;
    meanp.u += p.u;
  }
  meanp.e /= pts.size();
  meanp.i /= pts.size();
  meanp.u /= pts.size();
  // Find the stdev.
  for (const auto& p : pts) {
    stddevp.e += pow(p.e - meanp.e, 2.0);
    stddevp.i += pow(p.i - meanp.i, 2.0);
    stddevp.u += pow(p.u - meanp.u, 2.0);
  }
  stddevp.e = sqrt(stddevp.e / (pts.size() - 1));
  stddevp.i = sqrt(stddevp.i / (pts.size() - 1));
  stddevp.u = sqrt(stddevp.u / (pts.size() - 1));

  // printf("energy is %f,%f(%zu, mean %f s %f) and irmsd is %f,%f(%zu, mean %f s %f) and uirmsd is %f,%f(%zu, mean %f s %f)\n",
  //     minp.e, maxp.e, unique_e.size(), meanp.e, stddevp.e,
  //     minp.i, maxp.i, unique_i.size(), meanp.i, stddevp.i);
  //     minp.u, maxp.u, unique_u.size(), meanp.u, stddevp.u);

  // Star discrepancy, in all dims.
  float max_vol = vol(minp, maxp, stddevp);
  float max_disc = 0;
  point max_disc_p;
  int count = 0;

  for (const auto& p1 : pts) {
    // if (count++ % 100000 == 0) {
    //   printf("[%d/%d] disc %f\n", count, unique_e.size() * unique_i.size(), max_disc);
    // }
    float disc = one_disc(pts, meanp, p1, max_vol, stddevp);
    if (disc > max_disc) {
      max_disc = disc;
      max_disc_p = p1;
    }
  }
  printf("idx %d\n", max_disc_p.idx);
  return max_disc;
}

float centroid_discrepancy(const std::vector<point> &pts) {
  // Get the min/max in each direction.
  point minp = pts[0], maxp = pts[0], meanp, stddevp;
  for (const auto& p : pts) {
    minp.e = std::min(p.e, minp.e);
    minp.i = std::min(p.i, minp.i);
    minp.u = std::min(p.u, minp.u);
    maxp.e = std::max(p.e, maxp.e);
    maxp.i = std::max(p.i, maxp.i);
    maxp.u = std::max(p.u, maxp.u);

    meanp.e += p.e;
    meanp.i += p.i;
    meanp.u += p.u;
  }
  meanp.e /= pts.size();
  meanp.i /= pts.size();
  meanp.u /= pts.size();
  // Find the stdev.
  for (const auto& p : pts) {
    stddevp.e += pow(p.e - meanp.e, 2.0);
    stddevp.i += pow(p.i - meanp.i, 2.0);
    stddevp.u += pow(p.u - meanp.u, 2.0);
  }
  stddevp.e = sqrt(stddevp.e / (pts.size() - 1));
  stddevp.i = sqrt(stddevp.i / (pts.size() - 1));
  stddevp.u = sqrt(stddevp.u / (pts.size() - 1));

  // printf("energy is %f,%f(%zu, mean %f s %f) and irmsd is %f,%f(%zu, mean %f s %f) and uirmsd is %f,%f(%zu, mean %f s %f)\n",
  //     minp.e, maxp.e, unique_e.size(), meanp.e, stddevp.e,
  //     minp.i, maxp.i, unique_i.size(), meanp.i, stddevp.i);
  //     minp.u, maxp.u, unique_u.size(), meanp.u, stddevp.u);

  // Star discrepancy, in all dims.
  float max_vol = vol(minp, maxp, stddevp);
  float max_disc = 0;
  point max_disc_p;
  int count = 0;

  for (const auto& p1 : pts) {
    // if (count++ % 100000 == 0) {
    //   printf("[%d/%d] disc %f\n", count, unique_e.size() * unique_i.size(), max_disc);
    // }
    float disc = one_disc(pts, meanp, p1, max_vol, stddevp);
    if (disc > max_disc) {
      max_disc = disc;
      max_disc_p = p1;
    }
  }
  printf("idx %d\n", max_disc_p.idx);
  return max_disc;
}

float star_discrepancy(const std::vector<point> &pts) {
  // Get the min/max in each direction.
  point minp = pts[0], maxp = pts[0], meanp, stddevp;
  std::set<float> unique_e, unique_i, unique_u;
  for (const auto& p : pts) {
    minp.e = std::min(p.e, minp.e);
    minp.i = std::min(p.i, minp.i);
    minp.u = std::min(p.u, minp.u);
    maxp.e = std::max(p.e, maxp.e);
    maxp.i = std::max(p.i, maxp.i);
    maxp.u = std::max(p.u, maxp.u);

    unique_e.insert(p.e);
    unique_i.insert(p.i);
    unique_u.insert(p.u);

    meanp.e += p.e;
    meanp.i += p.i;
    meanp.u += p.u;
  }
  meanp.e /= pts.size();
  meanp.i /= pts.size();
  meanp.u /= pts.size();
  // Find the stdev.
  for (const auto& p : pts) {
    stddevp.e += pow(p.e - meanp.e, 2.0);
    stddevp.i += pow(p.i - meanp.i, 2.0);
    stddevp.u += pow(p.u - meanp.u, 2.0);
  }
  stddevp.e = sqrt(stddevp.e / (pts.size() - 1));
  stddevp.i = sqrt(stddevp.i / (pts.size() - 1));
  stddevp.u = sqrt(stddevp.u / (pts.size() - 1));

  // Need to get rid of the zero points as well.
  unique_e.erase(minp.e);
  unique_i.erase(minp.i);
  unique_u.erase(minp.u);

  // printf("energy is %f,%f(%zu, mean %f s %f) and irmsd is %f,%f(%zu, mean %f s %f) and uirmsd is %f,%f(%zu, mean %f s %f)\n",
  //     minp.e, maxp.e, unique_e.size(), meanp.e, stddevp.e,
  //     minp.i, maxp.i, unique_i.size(), meanp.i, stddevp.i);
  //     minp.u, maxp.u, unique_u.size(), meanp.u, stddevp.u);

  // Star discrepancy, in all dims.
  float max_vol = vol(minp, maxp, stddevp);
  float max_disc = 0;
  int count = 0;
  std::vector<float> evec(unique_e.begin(), unique_e.end());
#pragma omp parallel for reduction(max:max_disc)
  for (int ei = 0; ei < evec.size(); ++ei) {
    const auto& e = evec[ei];
    for (const auto& i : unique_i) {
      for (const auto& u : unique_u) {
        // if (count++ % 100000 == 0) {
        //   printf("[%d/%d] disc %f\n", count, unique_e.size() * unique_i.size(), max_disc);
        // }
        point test_p(e, i, u, -1);
        max_disc = std::max(max_disc, one_disc(pts, minp, test_p, max_vol, stddevp));
      }
    }
  }
  return max_disc;
}

int main(int argc, char* argv[]) {
  if (argc <= 5) {
    fprintf(stderr, "usage: discrepancy <input_file> <idx_col> <energy_col> <irmsd_col> <uirmsd_col>\n");
    return -1;
  }

  const char* inf = argv[1];
  int idxcol = atoi(argv[2]);
  int ecol = atoi(argv[3]);
  int icol = atoi(argv[4]);
  int ucol = atoi(argv[5]);

  auto points = read_inf(inf, idxcol, ecol, icol, ucol);

  // printf("Points have size %zu\n", points.size());
  float disc = centroid_discrepancy(points);
  printf("Discrepancy is %f\n", disc);
}
