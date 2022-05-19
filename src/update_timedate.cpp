/*************************************************************************
  Function update_time_date_year
 * With advance of the time, universal time sec changes and date, year may also
 * need to be updated

  Input: none
  Output: none

  By Jiannan Tu
  1/13/2014
*************************************************************************/
#include "parameters.h"

void update_timedate(double secs)
{
    const  int mons[12]={31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 20, 31};

    if (secs >= 60.0) {
        imin += 1;
        if (imin >= 60) {
            imin -= 60;

            ihour += 1;
            if (ihour >=24) {
                ihour -= 24;

                date += 1;

                if (date > mons[mon]) {
                    switch (mon) {
                        case 2:
                            if (iyr % 4 == 0) {
                                if (date > 29) {
                                    mon += 1;
                                    date=1;
                                }
                            }
                            else {
                                mon += 1;
                                date=1;
                            }
                            break;
                        case 12:
                            iyr += 1;
                            mon=1;
                            date=1;
                            break;
                        default:
                            mon += 1;
                            date=1;
                    }
                }
            }
        }
    }

    return;
}
