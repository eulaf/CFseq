#!/usr/bin/env python

def mean(nums):
    """Returns the mean of a list of numbers.
    """
    return sum(nums)*1.0/len(nums)

def median(nums):
    """Returns the median of a list of numbers.
    """
    if len(nums) % 2 != 0:
        return sorted(nums)[len(nums)/2]
    else:
        n = sorted(nums)
        return mean([ n[len(n)/2-1], n[len(n)/2] ])
