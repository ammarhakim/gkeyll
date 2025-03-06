#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define ANSI_COLOR_RED "\x1b[31m"
#define ANSI_COLOR_GREEN "\x1b[32m"
#define ANSI_COLOR_YELLOW "\x1b[33m"
#define ANSI_COLOR_BLUE "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN "\x1b[36m"
#define ANSI_COLOR_RESET "\x1b[0m"

#define RELATIVE_TOLERANCE pow(10.0, -16.0)

int system(const char *command);

void
runTest(const char* test_name, const char* test_name_human, const int test_output_count, const char test_outputs[][64])
{
  int counter = 0;

  char counter_buffer[128];
  snprintf(counter_buffer, 128, "ci/output/%s_counter.dat", test_name);
  FILE *counter_ptr = fopen(counter_buffer, "r");
  if (counter_ptr != NULL) {
    fscanf(counter_ptr, "%d", &counter);
    fclose(counter_ptr);
  }

  counter += 1;

  counter_ptr = fopen(counter_buffer, "w");
  fprintf(counter_ptr, "%d", counter);
  fclose(counter_ptr);

  printf("Running %s...\n", test_name_human);

  char command_buffer1[256];
  snprintf(command_buffer1, 256, "rm -rf ./%s-stat.json", test_name);
  system(command_buffer1);
  
  char command_buffer2[256];
  snprintf(command_buffer2, 256, "make build/regression/rt_%s > /dev/null 2>&1", test_name);
  system(command_buffer2);

  char command_buffer3[256];
  snprintf(command_buffer3, 256, "./build/regression/rt_%s -m > ./ci/output/rt_%s_%d.dat 2>&1", test_name, test_name, counter);
  system(command_buffer3);

  char file_buffer1[128];
  snprintf(file_buffer1, 128, "./%s-stat.json", test_name);
  FILE *file_ptr1 = fopen(file_buffer1, "r");
  if (file_ptr1 == NULL) {
    printf("*** Something catastrophic happened. Test aborting... ***\n");
  }
  else {
    char command_buffer4[256];
    snprintf(command_buffer4, 256, "mv ./%s-stat.json ci/output/%s-stat_%d.json", test_name, test_name, counter);
    system(command_buffer4);
  }

  for (int i = 0; i < test_output_count; i++) {
    char file_buffer2[128];
    snprintf(file_buffer2, 128, "./%s-%s.gkyl", test_name, test_outputs[i]);
    FILE *file_ptr2 = fopen(file_buffer2, "r");
    if (file_ptr2 == NULL) {
      printf("*** Something catastrophic happened. Test aborting... ***\n");
    }
    else {
      char command_buffer5[256];
      snprintf(command_buffer5, 256, "mv ./%s-%s.gkyl ci/output/%s-%s_%d.gkyl", test_name, test_outputs[i], test_name, test_outputs[i], counter);
      system(command_buffer5);
    }
  }

  printf("Finished %s.\n\n", test_name_human);
}

void
analyzeTestOutput(const char* test_name, const char* test_name_human, const int test_output_count, const char test_outputs[][64])
{
  printf("%s:\n\n", test_name_human);

  int counter = 0;

  char counter_buffer[64];
  snprintf(counter_buffer, 64, "ci/output/%s_counter.dat", test_name);
  FILE *counter_ptr = fopen(counter_buffer, "r");
  if (counter_ptr != NULL) {
    fscanf(counter_ptr, "%d", &counter);
    fclose(counter_ptr);
  }

  int failure[counter + 1];
  int updatecalls[counter + 1];
  int forwardeuler[counter + 1];
  int rk2failures[counter + 1];
  int rk3failures[counter + 1];
  double speciesrhs[counter + 1];
  double speciescollisionsrhs[counter + 1];
  double fieldrhs[counter + 1];
  double speciescollisionalmoments[counter + 1];
  double totalupdate[counter + 1];
  int memoryleakcount[counter + 1];
  char *memoryleaks[counter + 1];
  long double averages[counter + 1][test_output_count];

  for (int i = 1; i < counter + 1; i++) {
    char *output;
    long file_size;
    char buffer[128];
    snprintf(buffer, 128, "ci/output/rt_%s_%d.dat", test_name, i);

    FILE *output_ptr = fopen(buffer, "rb");
    fseek(output_ptr, 0, SEEK_END);
    file_size = ftell(output_ptr);
    rewind(output_ptr);
    output = calloc(file_size, (sizeof(char)));
    fread(output, sizeof(char), file_size, output_ptr);
    fclose(output_ptr);

    failure[i] = 0;

    updatecalls[i] = 0;
    if (strstr(output, "Number of update calls ") != NULL) {
      char *full_substring = strstr(output, "Number of update calls ");
      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Number of update calls ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Number of update calls ")];
        substring_index += 1;
      }

      char *end_ptr;
      updatecalls[i] = strtol(substring, &end_ptr, 10);
    }
    else {
      failure[i] = 1;
    }

    forwardeuler[i] = 0;
    if (strstr(output, "Number of forward-Euler calls ") != NULL) {
      char *full_substring = strstr(output, "Number of forward-Euler calls ");
      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Number of forward-Euler calls ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Number of forward-Euler calls ")];
        substring_index += 1;
      }

      char *end_ptr;
      forwardeuler[i] = strtol(substring, &end_ptr, 10);
    }
    else {
      failure[i] = 1;
    }

    rk2failures[i] = 0;
    if (strstr(output, "Number of RK stage-2 failures ") != NULL) {
      char *full_substring = strstr(output, "Number of RK stage-2 failures ");
      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Number of RK stage-2 failures ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Number of RK stage-2 failures ")];
        substring_index += 1;
      }

      char *end_ptr;
      rk2failures[i] = strtol(substring, &end_ptr, 10);
    }
    else {
      failure[i] = 1;
    }

    rk3failures[i] = 0;
    if (strstr(output, "Number of RK stage-3 failures ") != NULL) {
      char *full_substring = strstr(output, "Number of RK stage-3 failures ");
      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Number of RK stage-3 failures ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Number of RK stage-3 failures ")];
        substring_index += 1;
      }

      char *end_ptr;
      rk3failures[i] = strtol(substring, &end_ptr, 10);
    }
    else {
      failure[i] = 1;
    }

    speciesrhs[i] = 0.0;
    if (strstr(output, "Species RHS calc took ") != NULL) {
      char *full_substring = strstr(output, "Species RHS calc took ");
      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Species RHS calc took ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Species RHS calc took ")];
        substring_index += 1;
      }

      char *end_ptr;
      speciesrhs[i] = strtod(substring, &end_ptr);
    }
    else {
      failure[i] = 1;
    }

    speciescollisionsrhs[i] = 0.0;
    if (strstr(output, "Species collisions RHS calc took ") != NULL) {
      char *full_substring = strstr(output, "Species collisions RHS calc took ");
      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Species collisions RHS calc took ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Species collisions RHS calc took ")];
        substring_index += 1;
      }

      char *end_ptr;
      speciescollisionsrhs[i] = strtod(substring, &end_ptr);
    }
    else {
      failure[i] = 1;
    }

    fieldrhs[i] = 0.0;
    if (strstr(output, "Field RHS calc took ") != NULL) {
      char *full_substring = strstr(output, "Field RHS calc took ");
      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Field RHS calc took ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Field RHS calc took ")];
        substring_index += 1;
      }

      char *end_ptr;
      fieldrhs[i] = strtod(substring, &end_ptr);
    }
    else {
      failure[i] = 1;
    }

    speciescollisionalmoments[i] = 0.0;
    if (strstr(output, "Species collisional moments took ") != NULL) {
      char *full_substring = strstr(output, "Species collisional moments took ");
      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Species collisional moments took ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Species collisional moments took ")];
        substring_index += 1;
      }

      char *end_ptr;
      speciescollisionalmoments[i] = strtod(substring, &end_ptr);
    }
    else {
      failure[i] = 1;
    }

    totalupdate[i] = 0.0;
    if (strstr(output, "Total updates took ") != NULL) {
      char *full_substring = strstr(output, "Total updates took ");
      char substring[64];
      for (int j = 0; j < 64; j++) {
        substring[j] = '\0';
      }
      int substring_index = 0;

      while (full_substring[substring_index + strlen("Total updates took ")] != '\n') {
        substring[substring_index] = full_substring[substring_index + strlen("Total updates took ")];
        substring_index += 1;
      }

      char *end_ptr;
      totalupdate[i] = strtod(substring, &end_ptr);
    }
    else {
      failure[i] = 1;
    }
    
    if (failure[i] == 0) {
      char *temp = output;
      memoryleakcount[i] = 0;
      memoryleaks[i] = (char*)calloc(8192, sizeof(char));
      while (strstr(temp, "0x") != NULL) {
        temp = strstr(temp, "0x");

        char substring[64];
        for (int j = 0; j < 64; j++) {
          substring[j] = '\0';
        }

        int substring_index = 0;
        int valid_substring = 1;
        while (temp[substring_index] != ' ' && temp[substring_index] != '\n') {
          if (temp[substring_index] != '0' && temp[substring_index] != '1' && temp[substring_index] != '2' && temp[substring_index] != '3' && temp[substring_index] != '4'
            && temp[substring_index] != '5' && temp[substring_index] != '6' && temp[substring_index] != '7' && temp[substring_index] != '8' && temp[substring_index] != '9'
            && temp[substring_index] != 'a' && temp[substring_index] != 'b' && temp[substring_index] != 'c' && temp[substring_index] != 'd' && temp[substring_index] != 'e'
            && temp[substring_index] != 'f' && temp[substring_index] != 'x') {
            valid_substring = 0;
          }

          substring[substring_index] = temp[substring_index];
          substring_index += 1;
        }

        char *temp2 = output;
        int count = 0;
        while (strstr(temp2, substring) != NULL) {
          temp2 = strstr(temp2, substring);

          count += 1;
          temp2 += 1;
        }
        if (count == 1 && valid_substring == 1) {
          memoryleakcount[i] += 1;
          memoryleaks[i] = strcat(memoryleaks[i], substring);
          memoryleaks[i] = strcat(memoryleaks[i], " ");
        }
        
        temp += 1;
      }
    }

    for (int j = 0; j < test_output_count; j++) {
      char *data;
      long data_file_size;
      char data_buffer[256];
      snprintf(data_buffer, 256, "ci/output/%s-%s_%d.gkyl", test_name, test_outputs[j], i);

      FILE *data_ptr = fopen(data_buffer, "rb");

      if (data_ptr == NULL) {
        failure[i] = 1;
      }
      else {
        fseek(data_ptr, 0, SEEK_END);
        data_file_size = ftell(data_ptr);
        rewind(data_ptr);
        data = calloc(data_file_size, (sizeof(char)));
        fread(data, sizeof(char), data_file_size, data_ptr);
        fclose(data_ptr);

        long long total = 0;
        for (long k = 0; k < data_file_size; k++) {
          total += (long long)abs((int)data[k]);
        }

        averages[i][j] = (long double)total / (long double)data_file_size;
      }
    }
  }

  for (int i = 1; i < counter + 1; i++) {
    printf("Build number: %d\n", i);

    if (failure[i] == 1) {
      printf(ANSI_COLOR_RED "*** Catastrophic test failure ***" ANSI_COLOR_RESET "\n\n");
    }
    else {
      if (i == 1 || failure[i - 1] == 1) {
        printf("Update calls: %d\n", updatecalls[i]);
        printf("Forward-Euler calls: %d\n", forwardeuler[i]);
        printf("RK stage-2 failures: %d\n", rk2failures[i]);
        printf("RK stage-3 failures: %d\n", rk3failures[i]);
        printf("Species RHS time: %f\n", speciesrhs[i]);
        printf("Species collision RHS time: %f\n", speciescollisionsrhs[i]);
        printf("Field RHS time: %f\n", fieldrhs[i]);
        printf("Species collisional moments time: %f\n", speciescollisionalmoments[i]);
        printf("Total update time: %f\n", totalupdate[i]);
        if (memoryleakcount[i] != 0) {
          printf("Memory leaks: " ANSI_COLOR_RED "%s" ANSI_COLOR_RESET "\n", memoryleaks[i]);
        }
        else {
          printf("Memory leaks: " ANSI_COLOR_GREEN "None" ANSI_COLOR_RESET "\n");
        }
        printf("Correct: N/A\n\n");
      }
      else {
        if (updatecalls[i] != updatecalls[i - 1]) {
          printf("Update calls: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", updatecalls[i]);
        }
        else {
          printf("Update calls: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", updatecalls[i]);
        }

        if (forwardeuler[i] != forwardeuler[i - 1]) {
          printf("Forward-Euler calls: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", forwardeuler[i]);
        }
        else {
          printf("Forward-Euler calls: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", forwardeuler[i]);
        }

        if (rk2failures[i] > rk2failures[i - 1]) {
          printf("RK stage-2 failures: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", rk2failures[i]);
        }
        else {
          printf("RK stage-2 failures: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", rk2failures[i]);
        }

        if (rk3failures[i] > rk3failures[i - 1]) {
          printf("RK stage-3 failures: " ANSI_COLOR_RED "%d" ANSI_COLOR_RESET "\n", rk3failures[i]);
        }
        else {
          printf("RK stage-3 failures: " ANSI_COLOR_GREEN "%d" ANSI_COLOR_RESET "\n", rk3failures[i]);
        }

        if (speciesrhs[i] > speciesrhs[i - 1]) {
          if (speciesrhs[i - 1] > pow(10.0, -8.0)) {
            printf("Species RHS time: " ANSI_COLOR_RED "%f (+%.2f%%)" ANSI_COLOR_RESET "\n", speciesrhs[i],
              (((double)speciesrhs[i] / (double)speciesrhs[i - 1]) - 1.0) * 100.0);
          } else {
            printf("Species RHS time: " ANSI_COLOR_RED "%f (N/A)" ANSI_COLOR_RESET "\n", speciesrhs[i]);
          }
        }
        else {
          if (speciesrhs[i - 1] > pow(10.0, -8.0)) {
            printf("Species RHS time: " ANSI_COLOR_GREEN "%f (%.2f%%)" ANSI_COLOR_RESET "\n", speciesrhs[i],
              (((double)speciesrhs[i] / (double)speciesrhs[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Species RHS time: " ANSI_COLOR_GREEN "%f (N/A)" ANSI_COLOR_RESET "\n", speciesrhs[i]);
          }
        }

        if (speciescollisionsrhs[i] > speciescollisionsrhs[i - 1]) {
          if (speciescollisionsrhs[i - 1] > pow(10.0, -8.0)) {
            printf("Species collision RHS time: " ANSI_COLOR_RED "%f (+%.2f%%)" ANSI_COLOR_RESET "\n", speciescollisionsrhs[i],
              (((double)speciescollisionsrhs[i] / (double)speciescollisionsrhs[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Species collision RHS time: " ANSI_COLOR_RED "%f (N/A)" ANSI_COLOR_RESET, speciescollisionsrhs[i]);
          }
        }
        else {
          if (speciescollisionsrhs[i - 1] > pow(10.0, -8.0)) {
            printf("Species collision RHS time: " ANSI_COLOR_GREEN "%f (%.2f%%)" ANSI_COLOR_RESET "\n", speciescollisionsrhs[i],
              (((double)speciescollisionsrhs[i] / (double)speciescollisionsrhs[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Species collision RHS time: " ANSI_COLOR_GREEN "%f (N/A)" ANSI_COLOR_RESET "\n", speciescollisionsrhs[i]);
          }
        }

        if (fieldrhs[i] > fieldrhs[i - 1]) {
          if (fieldrhs[i - 1] > pow(10.0, -8.0)) {
            printf("Field RHS time: " ANSI_COLOR_RED "%f (+%.2f%%)" ANSI_COLOR_RESET "\n", fieldrhs[i],
              (((double)fieldrhs[i] / (double)fieldrhs[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Field RHS time: " ANSI_COLOR_RED "%f (N/A)" ANSI_COLOR_RESET "\n", fieldrhs[i]);
          }
        }
        else {
          if (fieldrhs[i - 1] > pow(10.0, -8.0)) {
            printf("Field RHS time: " ANSI_COLOR_GREEN "%f (%.2f%%)" ANSI_COLOR_RESET "\n", fieldrhs[i],
              (((double)fieldrhs[i] / (double)fieldrhs[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Field RHS time: " ANSI_COLOR_GREEN "%f (N/A)" ANSI_COLOR_RESET "\n", fieldrhs[i]);
          }
        }

        if (speciescollisionalmoments[i] > speciescollisionalmoments[i - 1]) {
          if (speciescollisionalmoments[i - 1] > pow(10.0, -8.0)) {
            printf("Species collisional moments time: " ANSI_COLOR_RED "%f (+%.2f%%)" ANSI_COLOR_RESET "\n", speciescollisionalmoments[i],
              (((double)speciescollisionalmoments[i] / (double)speciescollisionalmoments[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Species collisional moments time: " ANSI_COLOR_RED "%f (N/A)" ANSI_COLOR_RESET "\n", speciescollisionalmoments[i]);
          }
        }
        else {
          if (speciescollisionalmoments[i - 1] > pow(10.0, -8.0)) {
            printf("Species collisional moments time: " ANSI_COLOR_GREEN "%f (%.2f%%)" ANSI_COLOR_RESET "\n", speciescollisionalmoments[i],
              (((double)speciescollisionalmoments[i] / (double)speciescollisionalmoments[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Species collisional moments time: " ANSI_COLOR_GREEN "%f (N/A)" ANSI_COLOR_RESET "\n", speciescollisionalmoments[i]);
          }
        }

        if (totalupdate[i] > totalupdate[i - 1]) {
          if (totalupdate[i  - 1] > pow(10.0, -8.0)) {
            printf("Total update time: " ANSI_COLOR_RED "%f (+%.2f%%)" ANSI_COLOR_RESET "\n", totalupdate[i],
              (((double)totalupdate[i] / (double)totalupdate[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Total update time: " ANSI_COLOR_RED "%f (N/A)" ANSI_COLOR_RESET "\n", totalupdate[i]);
          }
        }
        else {
          if (totalupdate[i - 1] > pow(10.0, -8.0)) {
            printf("Total update time: " ANSI_COLOR_GREEN "%f (%.2f%%)" ANSI_COLOR_RESET "\n", totalupdate[i],
              (((double)totalupdate[i] / (double)totalupdate[i - 1]) - 1.0) * 100.0);
          }
          else {
            printf("Total update time: " ANSI_COLOR_GREEN "%f (N/A)" ANSI_COLOR_RESET "\n", totalupdate[i]);
          }
        }

        if (memoryleakcount[i] != 0) {
          printf("Memory leaks: " ANSI_COLOR_RED "%s" ANSI_COLOR_RESET "\n", memoryleaks[i]);
        }
        else {
          printf("Memory leaks: " ANSI_COLOR_GREEN "None" ANSI_COLOR_RESET "\n");
        }

        int correct = 1;
        for (int j = 0; j < test_output_count; j++) {
          if (fabsl(averages[i][j] - averages[i - 1][j]) > RELATIVE_TOLERANCE) {
            correct = 0;
          }
        }

        if ((updatecalls[i] != updatecalls[i - 1]) || (forwardeuler[i] != forwardeuler[i - 1]) ||
          (rk2failures[i] != rk2failures[i - 1]) || (rk3failures[i] != rk3failures[i - 1]) || (correct != 1)) {
          printf("Correct: " ANSI_COLOR_RED "No" ANSI_COLOR_RESET "\n\n");
        }
        else {
          printf("Correct: " ANSI_COLOR_GREEN "Yes" ANSI_COLOR_RESET "\n\n");
        }
      }
    }
  }
}

void
regenerateTest(const char* test_name, const int test_output_count, const char test_outputs[][64])
{
  int counter = 0;

  char counter_buffer[64];
  snprintf(counter_buffer, 64, "ci/output/%s_counter.dat", test_name);
  FILE *counter_ptr = fopen(counter_buffer, "r");
  if (counter_ptr != NULL) {
    fscanf(counter_ptr, "%d", &counter);
  }
  fclose(counter_ptr);

  for (int i = 1 ; i < counter + 1; i++) {
    char command_buffer[128];
    snprintf(command_buffer, 128, "rm -rf ci/output/rt_%s_%d.dat", test_name, i);
    system(command_buffer);

    char command_buffer2[128];
    snprintf(command_buffer2, 128, "rm -rf ci/output/%s-stat_%d.json", test_name, i);
    system(command_buffer2);

    for (int j = 0; j < test_output_count; j++) {
      char command_buffer3[256];
      snprintf(command_buffer3, 256, "rm -rf ci/output/%s-%s_%d.gkyl", test_name, test_outputs[j], i);
      system(command_buffer3);
    }
  }

  counter_ptr = fopen(counter_buffer, "w");
  fprintf(counter_ptr, "%d", 0);
  fclose(counter_ptr);
}

int
main(int argc, char **argv)
{
  // To add a new test:
  // 1. Increase test_count by 1.
  // 2. Place the name of the test in test_names (in alphabetical order),
  // bumping down other tests in groups of 5 (for ease of reading).
  // 3. Insert the human-readable form of the name in test_names_human.
  // 4. Insert the number of outputs to check in test_output_count.
  // 5. Insert the output file suffixes in test_outputs.

  const int test_count = 36;
  char test_names[test_count][64] = {
    // Please list in alphabetic order.
    "gk_asdex_3x2v_p1",
    "gk_bgk_3x2v_p1",
    "gk_bgk_relax_1x2v_p1",
    "gk_bgk_relax_bimaxwellian_1x2v_p1",
    "gk_bgk_relax_bimaxwellian_nonuniformv_1x2v_p1",
    //
    "gk_d3d_iwl_2x2v_p1",
    "gk_ion_sound_1x2v_p1",
    "gk_ion_sound_nonuniformv_1x2v_p1",
    "gk_ion_sound_adiabatic_elc_1x2v_p1",
    "gk_lapd_cart_3x2v_p1",
    //
    "gk_lapd_cyl_3x2v_p1",
    "gk_lbo_relax_1x2v_p1",
    "gk_lbo_relax_bimaxwellian_nonuniformv_1x2v_p1",
    "gk_lbo_relax_bimaxwellian_nonuniformv_3x2v_p1",
    "gk_lbo_relax_nonuniformv_1x2v_p1",
    //
    "gk_lbo_relax_varnu_1x2v_p1",
    "gk_li_react_3x2v_p1",
    "gk_ltx_1x2v_p1",
    "gk_ltx_boltz_elc_1x2v_p1",
    "gk_mdpx_cart_3x2v_p1",
    //
    "gk_mirror_boltz_elc_1x2v_p1",
    "gk_mirror_boltz_elc_nonuniformz_1x2v_p1",
    "gk_mirror_kinetic_elc_1x2v_p1",
    "gk_rad_1x2v_p1",
    "gk_rad_nonuniform_1x2v_p1",
    //
    "gk_sheath_1x2v_p1",
    "gk_sheath_2x2v_p1",
    "gk_sheath_3x2v_p1",
    "gk_sheath_bgk_1x2v_p1",
    "gk_sheath_nonuniformv_1x2v_p1",
    //
    "gk_sheath_nonuniformv_2x2v_p1",
    "gk_sheath_nonuniformv_3x2v_p1",
    "gk_solovev_out_3x2v_p1",
    "gk_step_out_2x2v_p1",
    "gk_wham_1x2v_p1",
    //
    "gk_wham_2x2v_p1",
  };
  char test_names_human[test_count][128] = {
    "3x2v_p1 ASDEX test",
    "3x2v_p1 BGK collision test",
    "1x2v_p1 BGK collision relaxation test",
    "1x2v_p1 BGK relaxation test w/ bimaxwellians",
    "1x2v_p1 BGK relaxation test w/ bimaxwellians and nonuniform v",
    //
    "2x2v_p1 DIII-D IWL test",
    "1x2v_p1 Ion sound test",
    "1x2v_p1 Ion sound test w/ nonuniform v",
    "1x2v_p1 Ion sound test w/ adiabatic electrons",
    "3x2v_p1 LAPD test in Cartesian coordinates",
    //
    "3x2v_p1 LAPD test in cylindrical coordinates",
    "1x2v_p1 LBO relaxation test",
    "1x2v_p1 LBO relaxation test w/ bimaxwellians and nonuniform v",
    "3x2v_p1 LBO relaxation test w/ bimaxwellians and nonuniform v",
    "1x2v_p1 LBO relaxation test w/ nonuniform v",
    //
    "1x2v_p1 LBO relaxation test w/ nu(x,t)",
    "3x2v_p1 Li reactions test",
    "1x2v_p1 LTX SOL test",
    "1x2v_p1 LTX SOL test w/ Boltzmann electrons",
    "3x2v_p1 MDPX test",
    //
    "1x2v_p1 mirror test w/ Boltzmann electrons",
    "1x2v_p1 mirror test w/ Boltzmann electrons and nonuniform z",
    "1x2v_p1 mirror test w/ kinetic electrons",
    "1x2v_p1 Radiation operator test",
    "1x2v_p1 Radiation operator test w/ nonuniform v",
    //
    "1x2v_p1 Sheath test",
    "2x2v_p1 Sheath test",
    "3x2v_p1 Sheath test",
    "1x2v_p1 Sheath test w/ BGK operator",
    "1x2v_p1 Sheath test w/ nonuniform v",
    //
    "2x2v_p1 Sheath test w/ nonuniform v",
    "3x2v_p1 Sheath test w/ nonuniform v",
    "3x2v_p1 Solovev (outboard SOL) test",
    "2x2v_p1 STEP (outboard SOL)test",
    "1x2v_p1 WHAM test",
    //
    "2x2v_p1 WHAM test",
  };
  int test_output_count[test_count] = { 
    14, 14, 8, 10, 14,
    14, 12, 16, 4, 14,
    14, 8, 10, 10, 12,
    8, 32, 14, 4, 23,
    8, 8, 14, 16, 20,
    14, 14, 14, 14, 18,
    18, 18, 14, 25, 14,
    14,
  };
  char test_outputs[test_count][64][64] = {
    {
      "elc_1", "elc_source_1", "elc_M0_1", "elc_M1_1", "elc_M2_1", "elc_integrated_moms",
      "ion_1", "ion_source_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "field_1", "field_energy"
    },
    {
      "elc_1", "elc_source_1", "elc_M0_1", "elc_M1_1", "elc_M2_1", "elc_integrated_moms",
      "ion_1", "ion_source_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "field_1", "field_energy"
    },
    { 
      "square_1", "square_M0_1", "square_M1_1", "square_integrated_moms",
      "bump_1", "bump_M0_1", "bump_M1_1", "bump_integrated_moms",
    },
    {
      "elc_1", "elc_M0_1", "elc_M1_1", "elc_M2par_1", "elc_M2perp_1",
      "ion_1", "ion_M0_1", "ion_M1_1", "ion_M2par_1", "ion_M2perp_1",
    },
    {
      "elc_1", "elc_M0_1", "elc_M1_1", "elc_M2par_1", "elc_M2perp_1", "elc_mapc2p_vel", "elc_jacobvel",
      "ion_1", "ion_M0_1", "ion_M1_1", "ion_M2par_1", "ion_M2perp_1", "ion_mapc2p_vel", "ion_jacobvel",
    },
//
    {
      "elc_1", "elc_source_1", "elc_M0_1", "elc_M1_1", "elc_M2_1", "elc_integrated_moms",
      "ion_1", "ion_source_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "field_1", "field_energy"
    },
    {
      "elc_1", "elc_M0_1", "elc_M1_1", "elc_M2_1", "elc_integrated_moms",
      "ion_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "field_1", "field_energy"
    },
    {
      "elc_1", "elc_M0_1", "elc_M1_1", "elc_M2_1", "elc_integrated_moms", "elc_mapc2p_vel", "elc_jacobvel",
      "ion_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "ion_mapc2p_vel", "ion_jacobvel", "field_1", "field_energy"
    },
    { 
      "ion_1", "ion_integrated_moms", "field_1", "field_energy"
    },
    {
      "elc_1", "elc_source_1", "elc_M0_1", "elc_M1_1", "elc_M2_1", "elc_integrated_moms",
      "ion_1", "ion_source_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "field_1", "field_energy",
    },
//
    {
      "elc_1", "elc_source_1", "elc_M0_1", "elc_M1_1", "elc_M2_1", "elc_integrated_moms",
      "ion_1", "ion_source_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "field_1", "field_energy",
    },
    { 
      "square_1", "square_M0_1", "square_M1_1", "square_integrated_moms",
      "bump_1", "bump_M0_1", "bump_M1_1", "bump_integrated_moms",
    },
    {
      "elc_1", "elc_M0_1", "elc_M1_1", "elc_M2par_1", "elc_M2perp_1",
      "ion_1", "ion_M0_1", "ion_M1_1", "ion_M2par_1", "ion_M2perp_1",
    },
    {
      "elc_1", "elc_M0_1", "elc_M1_1", "elc_M2par_1", "elc_M2perp_1",
      "ion_1", "ion_M0_1", "ion_M1_1", "ion_M2par_1", "ion_M2perp_1",
    },
    { 
      "square_1", "square_M0_1", "square_M1_1", "square_integrated_moms", "square_mapc2p_vel", "square_jacobvel",
      "bump_1", "bump_M0_1", "bump_M1_1", "bump_integrated_moms", "bump_mapc2p_vel", "bump_jacobvel",
    },
//
    { 
      "square_1", "square_M0_1", "square_M1_1", "square_integrated_moms",
      "bump_1", "bump_M0_1", "bump_M1_1", "bump_integrated_moms",
    },
    {
      "elc_1", "elc_source_1", "elc_M0_1", "elc_M1_1", "elc_M2par_1", "elc_M2perp_1", "elc_Li2_Li1_iz_react_1", "elc_Li2_Li1_recomb_react_1", "elc_integrated_moms",
      "ion_1", "ion_source_1", "ion_M0_1", "ion_M1_1", "ion_M2par_1", "ion_M2perp_1", "ion_integrated_moms", "field_1", "field_energy",
      "Li1_1", "Li1_source_1", "Li1_M0_1", "Li1_M1_1", "Li1_M2par_1", "Li1_M2perp_1", "Li1_integrated_moms",
      "Li2_1", "Li2_source_1", "Li2_M0_1", "Li2_M1_1", "Li2_M2par_1", "Li2_M2perp_1", "Li2_integrated_moms",
    },
    {
      "elc_1", "elc_source_1", "elc_M0_1", "elc_M1_1", "elc_M2_1", "elc_integrated_moms",
      "ion_1", "ion_source_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "field_1", "field_energy",
    },
    { 
      "ion_1", "ion_integrated_moms", "field_1", "field_energy"
    },
    {
      "elc_1", "elc_M0_1", "elc_M1_1", "elc_M2_1", "elc_integrated_moms", "elc_Ar1_Ar0_iz_react_neut_1", "elc_Ar1_Ar0_recomb_react_neut_1",
      "ion_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "field_1", "field_energy",
      "Ar0_0", "Ar0_M0_0", "Ar0_M1i_0", "Ar0_M2_0",
      "Ar1_1", "Ar1_M0_1", "Ar1_M1_1", "Ar1_M2_1", "Ar1_integrated_moms",
    },
//
    {
      "ion_1", "ion_source_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "field_1", "field_energy",
    },
    {
      "ion_1", "ion_source_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "field_1", "field_energy",
    },
    {
      "elc_1", "elc_source_1", "elc_M0_1", "elc_M1_1", "elc_M2_1", "elc_integrated_moms",
      "ion_1", "ion_source_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "field_1", "field_energy",
    },
    {
      "elc_1", "elc_source_1", "elc_M0_1", "elc_M1_1", "elc_M2_1", "elc_integrated_moms", "elc_adiation_integrated_moms",
      "ion_1", "ion_source_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "ion_emissivity_1", "field_1", "field_energy",
    },
    {
      "elc_1", "elc_source_1", "elc_M0_1", "elc_M1_1", "elc_M2_1", "elc_integrated_moms", "elc_mapc2p_vel", "elc_jacobvel", "elc_adiation_integrated_moms",
      "ion_1", "ion_source_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "ion_mapc2p_vel", "ion_jacobvel", "ion_emissivity_1", "field_1", "field_energy",
    },
//
    {
      "elc_1", "elc_source_1", "elc_M0_1", "elc_M1_1", "elc_M2_1", "elc_integrated_moms",
      "ion_1", "ion_source_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "field_1", "field_energy",
    },
    {
      "elc_1", "elc_source_1", "elc_M0_1", "elc_M1_1", "elc_M2_1", "elc_integrated_moms",
      "ion_1", "ion_source_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "field_1", "field_energy",
    },
    {
      "elc_1", "elc_source_1", "elc_M0_1", "elc_M1_1", "elc_M2_1", "elc_integrated_moms",
      "ion_1", "ion_source_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "field_1", "field_energy",
    },
    {
      "elc_1", "elc_source_1", "elc_M0_1", "elc_M1_1", "elc_M2_1", "elc_integrated_moms",
      "ion_1", "ion_source_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "field_1", "field_energy",
    },
    {
      "elc_1", "elc_source_1", "elc_M0_1", "elc_M1_1", "elc_M2_1", "elc_integrated_moms", "elc_mapc2p_vel", "elc_jacobvel",
      "ion_1", "ion_source_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "ion_mapc2p_vel", "ion_jacobvel", "field_1", "field_energy",
    },
//
    {
      "elc_1", "elc_source_1", "elc_M0_1", "elc_M1_1", "elc_M2_1", "elc_integrated_moms", "elc_mapc2p_vel", "elc_jacobvel",
      "ion_1", "ion_source_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "ion_mapc2p_vel", "ion_jacobvel", "field_1", "field_energy",
    },
    {
      "elc_1", "elc_source_1", "elc_M0_1", "elc_M1_1", "elc_M2_1", "elc_integrated_moms", "elc_mapc2p_vel", "elc_jacobvel",
      "ion_1", "ion_source_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "ion_mapc2p_vel", "ion_jacobvel", "field_1", "field_energy",
    },
    {
      "elc_1", "elc_source_1", "elc_M0_1", "elc_M1_1", "elc_M2_1", "elc_integrated_moms",
      "ion_1", "ion_source_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "field_1", "field_energy",
    },
    {
      "elc_1", "elc_M0_1", "elc_M1_1", "elc_M2_1", "elc_integrated_moms", "elc_Ar1_Ar0_iz_react_neut_1", "elc_Ar1_Ar0_recomb_react_neut_1", "elc_emissivity_Ar1_1", "elc_radiation_integrated_moms",
      "ion_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "field_1", "field_energy",
      "Ar0_0", "Ar0_M0_0", "Ar0_M1i_0", "Ar0_M2_0",
      "Ar1_1", "Ar1_M0_1", "Ar1_M1_1", "Ar1_M2_1", "Ar1_integrated_moms",
    },
    {
      "elc_1", "elc_source_1", "elc_M0_1", "elc_M1_1", "elc_M2_1", "elc_integrated_moms",
      "ion_1", "ion_source_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "field_1", "field_energy",
    },
//
    {
      "elc_1", "elc_source_1", "elc_M0_1", "elc_M1_1", "elc_M2_1", "elc_integrated_moms",
      "ion_1", "ion_source_1", "ion_M0_1", "ion_M1_1", "ion_M2_1", "ion_integrated_moms", "field_1", "field_energy",
    },
  };

  system("clear");
  system("mkdir -p ci/output");

  printf("** Gkeyll Gyrokinetics Automated Regression System **\n\n");

  if (argc > 1) {
    char *arg_ptr;

    if (strtol(argv[1], &arg_ptr, 10) == 1) {
      for (int i = 0; i < test_count; i++) {
        runTest(test_names[i], test_names_human[i], test_output_count[i], test_outputs[i]);
      }
    }
    else if (strtol(argv[1], &arg_ptr, 10) == 2) {
      for (int i = 0; i < test_count; i++) {
        analyzeTestOutput(test_names[i], test_names_human[i], test_output_count[i], test_outputs[i]);
      }
    }
    else if (strtol(argv[1], &arg_ptr, 10) == 3) {
      if (argc > 2) {
        if (strtol(argv[2], &arg_ptr, 10) >= 1 && strtol(argv[2], &arg_ptr, 10) <= test_count) {
          runTest(test_names[strtol(argv[2], &arg_ptr, 10) - 1], test_names_human[strtol(argv[2], &arg_ptr, 10) - 1],
            test_output_count[strtol(argv[2], &arg_ptr, 10) - 1], test_outputs[strtol(argv[2], &arg_ptr, 10) - 1]);
        }
        else {
          printf("Invalid test!\n");
        }
      }
      else {
        printf("Must specify which test to run!\n");
      }
    }
    else if (strtol(argv[1], &arg_ptr, 10) == 4) {
      if (argc > 2) {
        if (strtol(argv[2], &arg_ptr, 10) >= 1 && strtol(argv[2], &arg_ptr, 10) <= test_count) {
          analyzeTestOutput(test_names[strtol(argv[2], &arg_ptr, 10) - 1], test_names_human[strtol(argv[2], &arg_ptr, 10) - 1],
            test_output_count[strtol(argv[2], &arg_ptr, 10) - 1], test_outputs[strtol(argv[2], &arg_ptr, 10) - 1]);
        }
        else {
          printf("Invalid test!\n");
        }
      }
      else {
        printf("Must specify which test results to view!\n");
      }
    }
    else if (strtol(argv[1], &arg_ptr, 10) == 5) {
      for (int i = 0; i < test_count; i++) {
        regenerateTest(test_names[i], test_output_count[i], test_outputs[i]);
        runTest(test_names[i], test_names_human[i], test_output_count[i], test_outputs[i]);
      }
    }
    else if (strtol(argv[1], &arg_ptr, 10) == 6) {
      if (argc > 2) {
        if (strtol(argv[2], &arg_ptr, 10) >= 1 && strtol(argv[2], &arg_ptr, 10) <= test_count) {
          regenerateTest(test_names[strtol(argv[2], &arg_ptr, 10) - 1], test_output_count[strtol(argv[2], &arg_ptr, 10) - 1], test_outputs[strtol(argv[2], &arg_ptr, 10) - 1]);
          runTest(test_names[strtol(argv[2], &arg_ptr, 10) - 1], test_names_human[strtol(argv[2], &arg_ptr, 10) - 1],
            test_output_count[strtol(argv[2], &arg_ptr, 10) - 1], test_outputs[strtol(argv[2], &arg_ptr, 10) - 1]);
        }
        else {
          printf("Invalid test!\n");
        }
      }
      else {
        printf("Must specify which test results to (re)generate!\n");
      }
    }
    else {
      printf("Invalid option!\n");
    }
  }
  else {
    while (1) {
      printf("Please select an option to proceed:\n\n");
      printf("1 - Run Full Regression Suite\n");
      printf("2 - View All Regression Results\n");
      printf("3 - Run Specific Regression Test\n");
      printf("4 - View Specific Regression Result\n");
      printf("5 - (Re)generate All Accepted Results\n");
      printf("6 - (Re)generate Specific Accepted Result\n");
      printf("7 - Exit\n");

      int option;
      scanf("%d", &option);
      printf("\n");

      if (option == 1) {
        for (int i = 0; i < test_count; i++) {
          runTest(test_names[i], test_names_human[i], test_output_count[i], test_outputs[i]);
        }
      }
      else if (option == 2) {
        for (int i = 0; i < test_count; i++) {
          analyzeTestOutput(test_names[i], test_names_human[i], test_output_count[i], test_outputs[i]);
        }
      }
      else if (option == 3) {
        printf("Please select the test you wish to run:\n\n");
        for (int i = 0; i < test_count; i++) {
          printf("%d - %s\n", i + 1, test_names_human[i]);
        }

        int option2;
        scanf("%d", &option2);
        printf("\n");

        if (option2 >= 1 && option2 <= test_count) {
          runTest(test_names[option2 - 1], test_names_human[option2 - 1], test_output_count[option2 - 1], test_outputs[option2 - 1]);
        }
        else {
          printf("Invalid test!\n\n");
        }
      }
      else if (option == 4) {
        printf("Please select the test whose results you wish to view:\n\n");
        for (int i = 0; i < test_count; i++) {
          printf("%d - %s\n", i + 1, test_names_human[i]);
        }

        int option2;
        scanf("%d", &option2);
        printf("\n");

        if (option2 >= 1 && option2 <= test_count) {
          analyzeTestOutput(test_names[option2 - 1], test_names_human[option2 - 1], test_output_count[option2 - 1], test_outputs[option2 - 1]);
        }
        else {
          printf("Invalid test!\n\n");
        }
      }
      else if (option == 5) {
        for (int i = 0; i < test_count; i++) {
          regenerateTest(test_names[i], test_output_count[i], test_outputs[i]);
          runTest(test_names[i], test_names_human[i], test_output_count[i], test_outputs[i]);
        }
      }
      else if (option == 6) {
        printf("Please select the test whose accepted result you wish to (re)generate:\n\n");
        for (int i = 0; i < test_count; i++) {
          printf("%d - %s\n", i + 1, test_names_human[i]);
        }

        int option2;
        scanf("%d", &option2);
        printf("\n");

        if (option2 >= 1 && option2 <= test_count) {
          regenerateTest(test_names[option2 - 1], test_output_count[option2 - 1], test_outputs[option2 - 1]);
          runTest(test_names[option2 - 1], test_names_human[option2 - 1], test_output_count[option2 - 1], test_outputs[option2 - 1]);
        }
        else {
          printf("Invalid test!\n\n");
        }
      }
      else if (option == 7) {
        break;
      }
      else {
        printf("Invalid selection!\n\n");
      }
    }
  }

  return 0;
}
