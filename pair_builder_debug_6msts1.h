///////////////////////////////////////////////////////////////////////////////////////////////////////
//    This program receive hits from a febdriver and searches for ts1_ref events
//    Note: Set maximal timedifference of events to select wrt ts1_ref in the code
//    This Version is for permanent use (no counters)
//    This filter is for uBooNe!!!
//    Written by Thomas Mettler
///////////////////////////////////////////////////////////////////////////////////////////////////////

//#include <zmq.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <map>
//#include "uboone/CRT/CRTProducts/CRTHit.hh"
#include "CRTHit.hh"

#define MAX_TIME_PREBEAM 2000000 //ns before beam ts1
#define MAX_TIME_PASTBEAM 4000000 //ns after beam ts1
#define EVLEN 80        // event length of a raw event (80 for uBooNE)
#define WAIT 0       // wait x us after sending
#define EVSPERFEB 1024   // max events per feb per poll to buffer
#define MAXFEBNR 256
#define MSOVERLAP 10000000//50000
#define MAX_TIME_DIFFERENCE 400   //Set the maximal timedifference between hits

//define numbers to controll bufferstatus
#define PROBUF_READY_TO_FILL  0
#define PROBUF_FILLING        1
#define PROBUF_READY_TO_SCALE 2
#define PROBUF_SHIFT_SCALE    3
#define SCANBUF_READY_TO_FILL 0
#define SCANBUF_SCANNING      1

#define FILTER_PAIR_MODE 3
#define TS1_CORR 11

#define STRIPW 10.8

namespace crt{
  class pair_builder;
  struct EVENT_t;
  struct EVENT_t_send;
  struct EVENT_tpro;
  struct EOP_EVENT_t;
  struct SCAN_ORDER;
}

struct crt::EVENT_t{
		uint16_t mac5;
		uint16_t flags;
    uint16_t lostcpu;
    uint16_t lostfpga;
		uint32_t ts0;
		uint32_t ts1;
		uint16_t adc[32];
};
struct crt::EVENT_t_send {
		uint16_t mac5;
		uint16_t flags;
    uint16_t lostcpu;
    uint16_t lostfpga;
		uint32_t ts0;
		uint32_t ts1;
		uint16_t adc[32];
    uint16_t recover;
    uint32_t nrtrigger;
    uint32_t nrtrigger_11;
};
struct crt::EVENT_tpro {
		uint16_t mac5;
		uint16_t flags;
    uint16_t lostcpu;
		uint16_t lostfpga;
		uint32_t ts0;
		uint32_t ts1;
		uint16_t adc[32];
		uint32_t ts0_scaled;
    uint32_t ts1_scaled;
		uint32_t sec;
    uint32_t ts0_ref;
    uint32_t ms;
    uint16_t recover;
    uint32_t nrtrigger;
    uint32_t nrtrigger_11;
};
struct crt::EOP_EVENT_t {
		uint16_t mac5; // ==0xFFFF
		uint16_t flags; // ==0xFFFF
		uint16_t lostcpu;
		uint16_t lostfpga;
		uint32_t ts0; // ==MAGICWORD32
		uint32_t ts1; // ==MAGICWORD32
    int nevsinpoll; 
		uint32_t start_s;
		uint32_t d1;
		uint16_t start_ms;
		uint16_t dd2;
		uint32_t d2;
		uint32_t end_s;
		uint32_t d3;
		uint16_t end_ms;
};  // end-of-poll special event
struct crt::SCAN_ORDER{
  uint32_t sec;
  int ref_nr;
  uint32_t ts0_ref;
  int flags;
};

class crt::pair_builder{

public:
    
  //This will be our function call inside the module
  void make_pairs(const char * filename,const char * filename_store, int run_mode_, std::vector<crt::CRTHit>& allCRTHits);
  int find_pairs(const char * filename,const char * filename_store,int run_mode);
  void usage();   //gives you information how to run
  void receive_data(); // receive data from zmq socket
  void shift_scale(int mac, int ref_nr, int ts0_ref); //scale the timestamps and copy them for processing
  void scale_buffer();  // scale the timestams
  void scan_buffer_filter(int mac);  // scan/process all hits of the FEBs and searches for coincidences
  void filter_buffer(int mac);
  void scan_filter_buffer(int mac);
  unsigned int my_abs(unsigned int, unsigned int);  //define absolutevalue for uint
  
  
  int ev_counter_mac[MAXFEBNR+1];   //Number of events per module (mac) in the processing buffer
  int ev_counter_scan[MAXFEBNR+1];  //Number of events per module (mac) in the scanning buffer
  int ev_counter_filter_scan[MAXFEBNR+1];  //Number of events per module (mac) in the scanning buffer
  int ev_counter_filter[MAXFEBNR+1];
  
  uint32_t act_time[2][MAXFEBNR+1];    //number to read out the second and ms out of received special event [0]:sec, [1]:ms [][mac]:module [][MAXFEBNR]:time last poll
  uint32_t previous_sec, previous_ms;
  uint32_t previous2_sec;
  int event_time_diff[MAXFEBNR+1];
  int event_time_diff_old[MAXFEBNR+1];
  
  //Hit stuff
  int XY_pair(int order,int mac3, int mac4);
  int XYtracks_Top(int mac3, int mac4);
  int XYtracks_Pipe(int mac3, int mac4);
  int XYtracks_Feed(int mac3, int mac4);
  int XYtracks_Bottom(int mac3, int mac4);

  void InitReconst();
  double getFineCoordExperim(int sL, int sR);
  double TWCorrection(int max1ach, int max2ach);
  double getHitT1(int mac1, int strip1,  int t1, int mac2, int strip2, int t2);
  double getHitT2(int mac1, int strip1,  int t1, int mac2, int strip2, int t2);
  double getHitDT(int mac1, int strip1,  int t1, int max1_ach, int mac2, int strip2,int t2, int max2_ach);
  double getHitT(int mac1, int strip1,  int t1, int max1_ach, int mac2, int strip2,  int t2, int max2_ach);
  double getHitX(int mac1, int strip1, int mac2, int strip2);
  double getHitY(int mac1, int strip1, int mac2, int strip2);
  double getHitZ(int mac1, int strip1, int mac2, int strip2);
  
  void make2DHit(int);
  void check_storing(void);
  
  
  
  
  int ready_to_send[10], send_bufnr, ready_to_fill, ready_to_scan;
  FILE *data;
  long size_ev;
  long counter_tot;
  long counter_old;
  long stuck_event_counter;
  crt::SCAN_ORDER order_buffer[MAXFEBNR+1]; //to scan the events with the lowest second number
  //beamfilter variables
  int last2_ms[MAXFEBNR+1];
  int last1_ms[MAXFEBNR+1];
  int number_ms[MAXFEBNR+1];

  //beamfilter V2 variables
  int ts1ref_buffer[100][100];
  int ts1ref_counter[100];
  uint32_t ts1ref_second[100];
  uint32_t previous_sec_scan[MAXFEBNR+1];
  uint32_t previous_sec_mac[MAXFEBNR+1];
  uint32_t previous_sec_mac2[MAXFEBNR+1];

  int ts0ref_counter[MAXFEBNR+1];
  int nrwm1[MAXFEBNR+1], nrwm2[MAXFEBNR+1];
  int run_mode; //choose the mode (filtering, pairfinding etz...)

  int minus2_counter1, same_counter1;
  
  //crthit stuf
  int refused;
  double Xs[200][32];
  double Ys[200][32];
  double Zs[200][32];
  int Exists[200];
  double Dts[200];
  //const static double inch=2.54; //inch in cm

  int verbose_;

private:

  
  
};






