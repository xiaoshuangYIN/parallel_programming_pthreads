#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <string.h>
#include <sstream>
#include <utility>
#include <map>
#include <vector>
#include <pthread.h>
using namespace std;
int numThreads = 2;
struct timespec start_time, end_time;

typedef struct {
  int elev;
  vector< pair<int, int> > lowest_neibrs;
  int num_lowest_neibrs;
  double drop_abs;
  double drop_remained;
  double drop_to_trickle;
}point;

int  M;
int N;
double A;
double drop_sum;
point** matrix;
pthread_mutex_t ds_lock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t* arrayOfLocks;

double calc_time(struct timespec start, struct timespec end) {
  double start_sec = (double)start.tv_sec*1000000000.0 + (double)start.tv_nsec;
  double end_sec = (double)end.tv_sec*1000000000.0 + (double)end.tv_nsec;
  
  if (end_sec < start_sec) {
    return 0;
  } else {
    return end_sec - start_sec;
  }
}


/* function: find the lowest neighbor(s) for each point */
/* fill up the "lowest_neibrs" for each point in the matrix*/
void find_neighbors(int myId) {

  //compute bounds for this thread
  int startrow = myId * N/numThreads;
  int endrow = (myId+1) * (N/numThreads) - 1;
  
  //find neighbors over the strip of rows for this thread
  for (int i = startrow; i <= endrow; i++) {
    for (int j = 0; j < N; j++) {
      // initialize a vector to store all the  neighbors with elev and position
      vector<pair<int, pair<int, int> > > neibrs_vec;
      for (int x = -1; x < 2; x++) {
  	    for (int y = -1; y < 2; y++) {
  	      if ((x * y == 0) && // ignore nw, ne, sw, se
  	          (y + j >= 0) &&
  	          (x + i >= 0) &&
  	          (y + j < N) &&
  	          (x + i < N)) {
  	        pair<int, int> coordinates(i+x, j+y);
  	        pair<int, pair<int, int> > elev_coordinates(matrix[i+x][y+j].elev,coordinates);
  	        neibrs_vec.push_back(elev_coordinates);
  	      }
  	    }
      }
      //cout<<"length of vector: "<< neibrs_vec.size() << endl;
      int min = neibrs_vec[0].first;
      //find the lowest neighbor/neighbors
      for (unsigned int k = 0; k < neibrs_vec.size(); k++){
  	    //find the lower one
  	    if (neibrs_vec[k].first < min){
  	      min = neibrs_vec[k].first;
  	      matrix[i][j].lowest_neibrs.clear();
  	      matrix[i][j].lowest_neibrs.push_back(neibrs_vec[k].second);
  	    }
  	    //find the duplicate
  	    else if (neibrs_vec[k].first == min){
  	      matrix[i][j].lowest_neibrs.push_back(neibrs_vec[k].second);
  	    }
      }
      //the lowest position is the same of center point
      if (min == matrix[i][j].elev) {
  	    matrix[i][j].num_lowest_neibrs = 1;
  	    matrix[i][j].lowest_neibrs.clear();
  	    pair<int, int> itself(i,j);
  	    matrix[i][j].lowest_neibrs.push_back(itself);
      }
      // one or more neighbors is lower than the center point
      else {
  	    matrix[i][j].num_lowest_neibrs = matrix[i][j].lowest_neibrs.size();
      }
      //cout<<"num of lowest neighbirs : "<< matrix[i][j].lowest_neibrs.size()<<endl;
    }
  }

}

void *finder(void *arg) {
  int id = *((int*) arg);
  find_neighbors(id);
  free(arg);
  return NULL;
}

/* function: do the first traversal over all landscape points */
void *RecAbsTric(int myId) {

  
  //compute bounds for this thread
  int startrow = myId * N/numThreads;
  int endrow = (myId+1) * (N/numThreads) - 1;
  
  for (int i = startrow; i <= endrow; i++) {
   
    for (int j = 0; j < N; j++) {
	    // (1) receive one rain drop when it is raining
	    if (M > 0) {
	      matrix[i][j].drop_remained = matrix[i][j].drop_remained + 1;
	      //cout<< "remain = " << matrix[i][j].drop_remained <<endl;
	    }
	    //(2) absorb water when there is drop remained
	    if (matrix[i][j].drop_remained >= A){
	      matrix[i][j].drop_remained -= A;
	      matrix[i][j].drop_abs += A; 
        pthread_mutex_lock(&ds_lock);
	      drop_sum -= A; 
        pthread_mutex_unlock(&ds_lock);
	      //cout<< " drop_count = " << drop_sum <<endl;
	    }
	    else if (matrix[i][j].drop_remained > 0){
	      matrix[i][j].drop_abs += matrix[i][j].drop_remained; 
        pthread_mutex_lock(&ds_lock);
	      drop_sum -= matrix[i][j].drop_remained; 
        pthread_mutex_unlock(&ds_lock);
	      matrix[i][j].drop_remained  = 0;	
	    }
	    //(3) water to trickle
	    if (matrix[i][j].drop_remained >= 1){
	      matrix[i][j].drop_to_trickle = 1;
	      matrix[i][j].drop_remained -= 1;
	    }
	    else{
	      matrix[i][j].drop_to_trickle = matrix[i][j].drop_remained;
	      matrix[i][j].drop_remained  = 0;
	    }
    }
  }
  return NULL;
}

/* function: do the second traversal over all landscape points */
void *updateTrickle(int myId) {
  int startrow = myId * N/numThreads;
  int endrow = (myId+1) * (N/numThreads) - 1;
    
  for (int i = startrow; i <= endrow; i++) {
      for (int j = 0; j < N; j++) {
	    for (int k = 0; k < matrix[i][j].num_lowest_neibrs; k++){
          int coor_x = (matrix[i][j].lowest_neibrs)[k].first;
          int coor_y = (matrix[i][j].lowest_neibrs)[k].second;
          pthread_mutex_lock(&arrayOfLocks[coor_x * N + coor_y]);
	        matrix[coor_x][coor_y].drop_remained += 
	        (matrix[i][j].drop_to_trickle / matrix[i][j].num_lowest_neibrs);
          pthread_mutex_unlock(&arrayOfLocks[coor_x * N + coor_y]);
	    }
    }
  }
 
  return NULL;
}

void *rbt_worker(void *arg) {
  int id = *((int*) arg);
  RecAbsTric(id);
  free(arg);
  return NULL;
}

void *updater(void *arg) {
  int id = *((int*) arg);
  updateTrickle(id);
  free(arg);
  return NULL;
}

/* function: stimulate rain-fall run-off process */
int rain_fall() {
  drop_sum = 0.0;
  int time = 0;
  int *p;
  pthread_t *threads;
  int numOfLocks = N * N;
  arrayOfLocks = (pthread_mutex_t*) 
    malloc(numOfLocks * sizeof(pthread_mutex_t));
  for (int i=0; i < numOfLocks; i++) {
    arrayOfLocks[i] = PTHREAD_MUTEX_INITIALIZER;
  }
  
  //allocate thread handles
  threads = (pthread_t *) malloc(numThreads * sizeof(pthread_t));

  while(true){
    time++;
    printf("M = %d\n",M);
    printf("Drop sum = %0.3f\n",drop_sum);
    if (M > 0){
      drop_sum += N * N;
    }
    //create threads
    for (int i=0; i < numThreads; i++) {
      p = (int*) malloc(sizeof(int));
      *p = i;
      pthread_create(&threads[i], NULL, rbt_worker, (void*)(p));
    }
    //wait till all threads are done
    for (int i=0; i < numThreads; i++) {
      pthread_join(threads[i], NULL);
    }
    // stop condition!!!
    printf("%f\n", drop_sum);
    if (drop_sum <= 0 && M <= 0){
      free(threads);
      return time;
    }
    M--;
    //tranverse again to update neighbors
    for (int i=0; i < numThreads; i++) {
      p = (int*) malloc(sizeof(int));
      *p = i;
      pthread_create(&threads[i], NULL, updater, (void*)(p));
    }
    for (int i=0; i < numThreads; i++) {
      pthread_join(threads[i], NULL);
    }
  }
}

/* function: print out the elev of each node in the matrx */
void print_matrix_elev(point** matrix) {
  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++){
      cout << matrix[i][j].elev << " ";
    }
    cout << endl;
  }
}
    
int main(int argc, char** argv){  
  if (argc != 5) {
    cerr << "Usage: ./rainfall <M> <A> <N> <elevation_file>\n";
    return EXIT_FAILURE;
  }
  
  clock_gettime(CLOCK_MONOTONIC, &start_time);  


  //read in arguments
  ifstream elev_file (argv[4]);
  M = stoi(argv[1]);
  A = strtof(argv[2], NULL);
  N = stoi(argv[3]);
  //initialize matrix
  matrix = (point**) malloc(N * sizeof(point*));
  for (int i = 0; i < N; i++) {
    matrix[i] = (point*) malloc(N * sizeof(point));
    for (int j = 0; j < N; j++) {
      matrix[i][j].elev = 0;
      matrix[i][j].drop_remained = 0;
      matrix[i][j].drop_abs = 0.0;
      matrix[i][j].drop_to_trickle = 0.0;
      //matrix[i][j].lowset_neibrs = NULL;
      matrix[i][j].num_lowest_neibrs = 0;
    }
  }
  //check open file succeed or not
  if (!elev_file.is_open()) {
    cerr <<"CAN NOT OPEN FILE" <<endl;
    return EXIT_FAILURE;
  }
  //read file
  int num_row = 0;
  int num_col = 0;
  string line;

  while (getline(elev_file, line)){
    if (line.empty()){
      continue;
    }
    
    num_col = 0;
    string buf;
    stringstream stream(line);
    while (stream >> buf){
      if (num_col >= N) {
        cout<<"number of columns exceeds N\n";
        return EXIT_FAILURE;
      }
      matrix[num_row][num_col].elev = stoi(buf);
      num_col++;
    }

    if (num_col != N){
      cout<<"Row "<<num_row
          << " has " << num_col
          << "numbers, it is suppose to be "
          << N <<endl;
      return EXIT_FAILURE;
    }
    num_row += 1;
  }

  if (num_row != N){
    cout << "The number of rows of input matrix is not N\n";
    return EXIT_FAILURE;
  }

  elev_file.close();
  //print_matrix_elev(matrix);
  //cout<<"read matrix finished\n";

  //Pthread
  pthread_t *threads;
  int *p;
  
  //allocate thread handles
  threads = (pthread_t*) malloc(numThreads * sizeof(pthread_t));
  
  //create threads
  for (int i=0; i < numThreads; i++) {
    p = (int*) malloc(sizeof(int));
    *p = i;
    pthread_create(&threads[i], NULL, finder, (void*)(p));
  }
  //wait till all threads are done
  for (int i=0; i < numThreads; i++) {
    pthread_join(threads[i], NULL);
  }
  
  int time = rain_fall();
  

  cout << "Rainfall simulation took "<< time << " time steps to complete." << endl;

  clock_gettime(CLOCK_MONOTONIC, &end_time);
  // print time
  double elapsed_ns = calc_time(start_time, end_time);
  double elapsed_ms = elapsed_ns / 1000000.0;
  double elapsed_s =  elapsed_ms / 1000.0;
  printf("Time=%f milliseconds\n", elapsed_ms);
  printf("Time=%f seconds\n", elapsed_s);
  
   
  cout << "The following grid shows the number of raindrops absorbed at each point:" << endl;
  unsigned unsignedN = (unsigned int) N;
  for (unsigned i=0; i < unsignedN; i++) {
      for (unsigned j=0; j < unsignedN; j++) {
        cout << setw(8) << setprecision(6) << matrix[i][j].drop_abs;
      } //for j
      
      cout << endl;
   } //for i
  
    
  //free malloc's
  free(arrayOfLocks);
  free(threads);
  for (int i = 0; i < N; i++) {
    free(matrix[i]);
  }
  free(matrix);

  
 
  return 0;
}
