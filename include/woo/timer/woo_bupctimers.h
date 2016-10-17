/**
 *  Project:
 *
 *  File: woo_bupctimers.h
 *  Created: Feb 23, 2014
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <upc_tick.h>

typedef struct {
	upc_tick_t start_;
	upc_tick_t stop_;
	uint64_t elapsed_;
	int is_running_;
	int is_paused_;
} UPCTimer;


void upctimer_reset(UPCTimer* t) {
	t->start_ = 0;
	t->stop_ = 0;
	t->elapsed_ = 0;
	t->is_running_ = 0;
	t->is_paused_ = 0;
} // upctimer_reset()


void upctimer_start(UPCTimer* t) {
	upctimer_reset(t);
	t->is_running_ = 1;
	t->start_ = upc_ticks_now();
} // upctimer_start()


void upctimer_stop(UPCTimer* t) {
	if(t->is_running_ == 0) {
		fprintf(stderr, "error: timer not running\n");
		return;
	} // if
	if(t->is_paused_ == 0) t->stop_ = upc_ticks_now();
	t->is_running_ = 0;
	t->elapsed_ = upc_ticks_to_ns(t->stop_ - t->start_);
} // upctimer_stop()


void upctimer_pause(UPCTimer* t) {
	if(t->is_running_ == 0) {
		fprintf(stderr, "error: timer not running\n");
		return;
	} // if
	if(t->is_paused_ == 1) {
		fprintf(stderr, "warning: timer is already paused. ignoring\n");
		return;
	} // if
	t->stop_ = upc_ticks_now();
	t->is_paused_ = 1;
} // upctimer_pause()


void upctimer_resume(UPCTimer* t) {
	if(t->is_running_ == 0) {
		fprintf(stderr, "error: timer not running\n");
		return;
	} // if
	if(t->is_paused_ == 0) {
		fprintf(stderr, "warning: timer is not paused. ignoring\n");
		return;
	} // if
	upc_tick_t temp1 = upc_ticks_now();
	upc_tick_t temp2 = temp1 - t->stop_;
	t->start_ = t->start_ + temp2;
	t->is_paused_ = 0;
} // upctimer_resume()


double upctimer_elapsed_nsec(UPCTimer* t) {
	return (double) t->elapsed_;
} // upctimer_elapsed_nsec()


double upctimer_elapsed_usec(UPCTimer* t) {
	return (double) t->elapsed_ / 1e3;
} // upctimer_elapsed_nsec()


double upctimer_elapsed_msec(UPCTimer* t) {
	return (double) t->elapsed_ / 1e6;
} // upctimer_elapsed_nsec()


double upctimer_elapsed_sec(UPCTimer* t) {
	return (double) t->elapsed_ / 1e9;
} // upctimer_elapsed_nsec()
