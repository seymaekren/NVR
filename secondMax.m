function retval = secondMax(v)

[values] = sort(-v);
retval = -values(2);