#pragma once

void* ermod_connect();
void ermod_disconnect(void* handle);
int ermod_send_trajctory(void* handle_p, int natoms, double cell[], double coords[]);
