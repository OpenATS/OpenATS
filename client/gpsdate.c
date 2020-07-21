/*
 * Simple tool to convert GPS time and calendar dates.
 */

#include <stdio.h>
#include <stdlib.h>

/*
 * Return Modified Julian Day given calendar year,
 * month (1-12), and day (1-31).
 * - Valid for Gregorian dates from 17-Nov-1858.
 * - Adapted from sci.astro FAQ.
 */

long
DateToMjd (long Year, long Month, long Day)
{
    return
        367 * Year
        - 7 * (Year + (Month + 9) / 12) / 4
        - 3 * ((Year + (Month - 9) / 7) / 100 + 1) / 4
        + 275 * Month / 9
        + Day
        + 1721028
        - 2400000;
}

/*
 * Convert Modified Julian Day to calendar date.
 * - Assumes Gregorian calendar.
 * - Adapted from Fliegel/van Flandern ACM 11/#10 p 657 Oct 1968.
 */

void
MjdToDate (long Mjd, long *Year, long *Month, long *Day)
{
    long J, C, Y, M;

    J = Mjd + 2400001 + 68569;
    C = 4 * J / 146097;
    J = J - (146097 * C + 3) / 4;
    Y = 4000 * (J + 1) / 1461001;
    J = J - 1461 * Y / 4 + 31;
    M = 80 * J / 2447;
    *Day = J - 2447 * M / 80;
    J = M / 11;
    *Month = M + 2 - (12 * J);
    *Year = 100 * (C - 49) + Y + J;
}

/*
 * Convert GPS Week and Seconds to Modified Julian Day.
 * - Ignores UTC leap seconds.
 */

long GpsToMjd (long GpsCycle, long GpsWeek, long GpsSeconds)
{
    long GpsDays;

    GpsDays = ((GpsCycle * 1024) + GpsWeek) * 7 + (GpsSeconds / 86400);
    return DateToMjd(1980, 1, 6) + GpsDays;
}

/*
 * Test program.
 */

main (int argc, char *argv[])
{
    long Year, Month, Day;
    long Mjd;
    long GpsCycle, GpsWeek, GpsSeconds;

    switch (argc - 1) {
    /* Given MJD */
    case 1 :
        Mjd = atol(argv[1]);

        MjdToDate(Mjd, &Year, &Month, &Day);

        printf("Mjd %ld = %.4ld-%.2ld-%.2ld\n", Mjd, Year, Month, Day);
        break;

    /* Given GPS week and seconds */
    case 2 :
        GpsWeek = atol(argv[1]);
        GpsSeconds = atol(argv[2]);
        for (GpsCycle = 0; GpsCycle < 3; GpsCycle += 1) {

            Mjd = GpsToMjd(GpsCycle, GpsWeek, GpsSeconds);
            MjdToDate(Mjd, &Year, &Month, &Day);

            printf("Gps Cycle %ld, Week %ld, Seconds %ld ",
                   GpsCycle, GpsWeek, GpsSeconds);
            printf(" = Mjd %ld = %.4ld-%.2ld-%.2ld\n",
                   Mjd, Year, Month, Day);
        }
        break;

    /* Given Calendar year, month, and day */
    case 3:
        Year = atol(argv[1]);
        Month = atol(argv[2]);
        Day = atol(argv[3]);

        Mjd = DateToMjd (Year, Month, Day);

        printf("%.4ld-%.2ld-%.2ld = Mjd %ld\n", Year, Month, Day, Mjd);
        break;

    default :
        fprintf(stderr, "\nThis tool converts calendar date, MJD, and GPS dates.\n\n");
        fprintf(stderr, "Use %s [Mjd]\n", argv[0]);
        fprintf(stderr, "Use %s [Gps Week] [Seconds]\n", argv[0]);
        fprintf(stderr, "Use %s [Year] [Month] [Day]\n", argv[0]);
        break;
    }
}
