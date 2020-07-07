/*
 * Copyright 2020 Xiongjie Yu, Ryan Levy, and Bryan K. Clark
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once
#include <unistd.h> //for getpid()
int make_new_seed()
{

  time_t tp;
  struct tm *temp;
  unsigned int temp2, temp3;
  static unsigned int temp4 = 0xe0e1;

  time(&tp);
  temp = localtime(&tp);

  temp2 = (temp->tm_sec<<26)+(temp->tm_min<<20)+(temp->tm_hour<<15)+
    (temp->tm_mday<<10)+(temp->tm_mon<<6);
  temp3 = (temp->tm_year<<13)+(temp->tm_wday<<10)+(temp->tm_yday<<1)+
    temp->tm_isdst;
  temp2 ^= clock()^temp3;

  temp4 = (temp4*0xeeee)%0xffff;
  temp2 ^= temp4<<16;
  temp4 = (temp4*0xaeee)%0xffff;
  temp2 ^= temp4;

  temp2 &= 0x7fffffff;

  return temp2;
}

long make_seed(){
  time_t seconds;
  long s, seed, pid;

  pid = getpid();
  s = time ( &seconds );

  seed = abs(((s*181)*((pid-83)*359))%104729);
  return seed;
}
