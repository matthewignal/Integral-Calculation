
/*
 * Given the actual number of the function calls, count, and the total cpu 
 * time spent, time, estimate the number of function calls required so 
 * that the total cpu time consumed by the function, is between tmin and tmax. 
 * All times - time, tmin, and tmax - must be in the same units, e.g. secs,
 * usecs, mins, etc)
 */
int adjust_rep_count (int count, double time, double tmin, double tmax);

int adjust_rep_count (int count, double time, double tmin, double tmax)
{
    if (time > tmax)
    {
        count /= time/tmax;
        count = count < 1 ? 1 : count;
    }
    else if (time < tmin)
    {
        count /= time/tmin;
    }
    return count;
}
