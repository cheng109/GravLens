//
//  main.cpp
//  junLensCode
//
//  Created by cheng109 on 9/5/15.
//  Copyright (c) 2015 cheng109. All rights reserved.
//

#include <stdio.h>
//#include "image.h"
#include "fitsmanipulation.h"

int main(int argc, const char * argv[]) {
    // insert code here...
    char filename[]="image.fits.gz";
    readFits(filename);
    printf("hello world, \n");
    return 0;

}




