#include <algorithm>
#include <chrono>       // std::chrono::system_clock
#include <cstdio>
#include <cstring>
#include <random>       // std::default_random_engine
#include <vector>

#define DIFFS_6D 0.9689

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

std::vector<float> read_inf(const char* inf, int idxcol, int colno) {
  std::vector<float> ret;
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
    float val;
    if (!get_col(line, colno, &val)) {
      fprintf(stderr, "Error on line %d: Could not find value column %d!! on line %s\n", line_num, colno, line);
    }
    // Ignore the input.
    if (idx == 9999) continue;

    ret.push_back(val);
  }

  fclose(fp);
  if (line) free(line);

  return ret;
}

// vals is not const so we can shuffle it.
std::vector<std::vector<float>> prob_diff_chernoff(
    const std::vector<float>& eval_points, std::vector<float>& vals,
    float global_min, int point_size, int num_draws, unsigned seed) {
  std::vector<std::vector<float>> bounds_per_run;
  for (int run = 0; run < num_draws; ++run) {
    std::vector<float> bounds(eval_points.size(), 0);
    std::shuffle(vals.begin(), vals.end(), std::default_random_engine(seed));

    // Get `point_size` things and compute the Chernoff bounds.
    for (int bi = 0; bi < eval_points.size(); ++bi) {
      for (int i = 0; i < point_size; ++i) {
        if (vals[i] < eval_points[bi]) {
          bounds[bi]++;
        }
      }
      // Then normalize
      bounds[bi] /= point_size;
    }
    /*
    for (float b : bounds) {
      printf("%f ", b);
    }
    printf("\n");
    */
    //printf("%d %f - %f vs %f\n", point_size, min_run-global_min, min_run, global_min);
    //printf("%d %f %f\n", point_size, min_run, min_run-global_min);
    bounds_per_run.push_back(bounds);
  }
  
  return bounds_per_run;
}
  
float l2norm(const std::vector<float>& a, const std::vector<float>& b) {
  float total = 0;
  for (int i = 0; i < a.size(); ++i) {
    total += pow(a[i] - b[i], 2.0);
  }
  return pow(total, 0.5);
}

void run_chernoff_bounds(
    std::vector<float>& vals, float global_min,
    unsigned seed, int num_reps, int stepsize) {
  const std::vector<float> eval_points = {10, 9, 8, 7, 6, 5, 4, 3, 2, 1};
  auto run_results_full_data_ex = prob_diff_chernoff(
      eval_points, vals, global_min, vals.size(), 1, seed);
  std::vector<float> run_results_full_data = run_results_full_data_ex[0];

  // Initial with #stepsize points.
  auto run_results_full = prob_diff_chernoff(
      eval_points, vals, global_min, stepsize, 1 /* num_draws */, seed);
  std::vector<float> run_results_prev = run_results_full[0];
  float l2_ex_init = l2norm(run_results_prev, run_results_full_data);
  printf("%d 1.0 %f\n", stepsize, l2_ex_init);

  for (int point_size = stepsize*2; point_size <= vals.size(); point_size += stepsize) {
    run_results_full = prob_diff_chernoff(eval_points, vals, global_min, point_size, 1, seed);
    std::vector<float> run_results = run_results_full[0];

    float l2 = l2norm(run_results, run_results_prev);
    float l2_ex = l2norm(run_results, run_results_full_data);
    printf("%d %f %f\n", point_size, l2, l2_ex);
    run_results_prev = run_results;
  }
}

// vals is not const so we can shuffle it.
float prob_diff(
    std::vector<float>& vals, float global_min,
    int point_size, int num_draws, unsigned seed) {
  std::vector<float> mins_per_run;
  for (int run = 0; run < num_draws; ++run) {
    std::shuffle(vals.begin(), vals.end(), std::default_random_engine(seed));

    // Get `point_size` things and see how far the min is from the global min.
    float min_run = vals[0];
    for (int i = 1; i < point_size; ++i) {
      min_run = std::min(min_run, vals[i]);
    }
    mins_per_run.push_back(min_run);
    //printf("%d %f - %f vs %f\n", point_size, min_run-global_min, min_run, global_min);
    printf("%d %f %f\n", point_size, min_run, min_run-global_min);
  }
  
  float avg_val = 0;
  for (int i = 0; i < num_draws; ++i) {
    avg_val += (mins_per_run[i] - global_min);
  }
  return avg_val / num_draws;
}

void run_discrepancy_runs(std::vector<float>& vals, float global_min, unsigned seed, int num_reps, int stepsize) {
  float prev_diff = prob_diff(vals, global_min, stepsize, num_reps, seed);
  int num_same = 0;
  for (int i = stepsize*2; i <= vals.size(); i+=stepsize) {
    //printf("%d %f", i, prob_diff(vals, 
    float diff = prob_diff(vals, global_min, i, num_reps, seed);
    if (std::abs(prev_diff - diff) < 0.01) {
      fprintf(stderr, "[%d] %d %f-%f=%f\n", num_same, i, prev_diff, diff, std::abs(prev_diff-diff));
      if (++num_same > 5) {
        break;
      }
    } else {
      num_same = 0;
    }
    prev_diff = diff;
  }
}

int main(int argc, char* argv[]) {
  if (argc < 4) {
    fprintf(stderr, "usage: prob_diff <input> <index_column_no> <column_no> [step_size=10]\n");
    return -1;
  }
  char* inf = argv[1];
  int idxcol = atoi(argv[2]);
  int colno = atoi(argv[3]);
  int stepsize = 10;
  if (argc > 4) {
    stepsize = atoi(argv[4]);
  }

  int num_reps=20;

  // obtain a time-based seed:
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

  std::vector<float> vals = read_inf(inf, idxcol, colno);

  // Get the min element.
  float global_min = *std::min_element(vals.begin(), vals.end());

  //run_discrepancy_runs(vals, global_min, seed, num_reps, stepsize);
  run_chernoff_bounds(vals, global_min, seed, num_reps, stepsize);
}

