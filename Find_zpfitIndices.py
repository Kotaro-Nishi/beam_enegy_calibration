import matplotlib.pyplot as plt

import UndulatorDiffractionOneDimensional_more_U1U2 as UnDi
import InterferenceFunction_more as IFunc
import AmplitudeFormula as AUndu

from numpy import transpose as nptp
import numpy as np

import cv2
import os , sys
import re

import csv

from natsort import natsorted
from scipy.optimize import curve_fit
from lmfit import minimize, Parameters, report_fit, Minimizer

import random

LengthOfInterval = 425
NumOfSeries = 100
random.seed(10000)

#0,1,2,3,4,5
List000 = [5, 7, 11, 18, 20, 22, 25, 28, 29, 37, 38, 41, 48, 55, 65, 67, 68, 72, 75, 76, 78, 79, 82, 84, 85, 96, 100, 103, 104, 106, 109, 114, 119, 129, 131, 133, 139, 142, 146, 147, 155, 164, 165, 166, 167, 173, 175, 180, 181, 194, 202, 213, 214, 215, 217, 220, 221, 225, 228, 229, 232, 237, 241, 242, 244, 250, 252, 255, 263, 265, 267, 274, 279, 280, 282, 286, 289, 290, 291, 294, 296, 299, 301, 302, 303, 308, 309, 310, 312, 315, 328, 335, 338, 344, 346, 348, 355, 356, 358, 359, 367, 373, 374, 377, 380, 382, 383, 387, 388, 392, 407, 409, 423, 430, 435, 436, 438, 442, 443, 450, 451, 452, 454, 455, 458, 460, 463, 466, 468, 479, 482, 483, 489, 495, 499]

i_folder = 0

listOfboots = []

for ii in range(500):
    if ii not in List000:
        for i in range(NumOfSeries):
            random.randint(0,LengthOfInterval-1)
    else:
        for i in range(NumOfSeries):
            listOfboots.append(random.randint(i_folder*25,i_folder*25+LengthOfInterval-1))
            
num_bins = 425
n, bins, patches = plt.hist(listOfboots, num_bins, facecolor = 'green')

print(bins)
plt.show()

List025 = [2, 3, 8, 9, 10, 11, 13, 14, 15, 17, 19, 21, 22, 23, 24, 27, 34, 35, 37, 38, 39, 42, 43, 44, 46, 48, 49, 50, 52, 53, 54, 57, 58, 63, 64, 71, 73, 76, 78, 79, 81, 82, 85, 88, 93, 94, 95, 96, 97, 98, 101, 106, 107, 115, 116, 118, 120, 121, 122, 124, 125, 126, 127, 132, 134, 139, 141, 143, 147, 149, 150, 151, 154, 155, 156, 160, 162, 163, 165, 166, 168, 169, 171, 172, 176, 178, 189, 191, 193, 195, 196, 197, 198, 200, 203, 206, 208, 209, 213, 214, 215, 221, 224, 225, 227, 228, 229, 230, 232, 233, 234, 237, 239, 240, 241, 242, 243, 246, 247, 248, 249, 251, 256, 257, 258, 259, 261, 262, 265, 267, 268, 270, 271, 273, 276, 277, 278, 279, 280, 281, 282, 283, 285, 286, 287, 288, 289, 290, 292, 293, 297, 298, 302, 304, 305, 306, 307, 309, 310, 311, 313, 315, 318, 321, 322, 324, 326, 329, 330, 331, 332, 335, 337, 338, 339, 341, 346, 348, 350, 353, 354, 356, 357, 359, 360, 362, 364, 366, 367, 369, 370, 373, 376, 377, 378, 379, 382, 385, 387, 391, 392, 393, 394, 396, 397, 399, 401, 402, 403, 404, 408, 409, 410, 411, 413, 414, 416, 417, 419, 420, 421, 423, 424, 425, 427, 428, 429, 435, 436, 438, 439, 440, 442, 443, 446, 448, 449, 455, 457, 458, 461, 463, 465, 470, 471, 473, 474, 478, 482, 483, 491, 492, 493]

i_folder = 1

listOfboots = []

for ii in range(500):
    if ii not in List025:
        for i in range(NumOfSeries):
            random.randint(0,LengthOfInterval-1)
    else:
        for i in range(NumOfSeries):
            listOfboots.append(random.randint(i_folder*25,i_folder*25+LengthOfInterval-1))
            
n, bins, patches = plt.hist(listOfboots, num_bins, facecolor = 'green')
print(bins)
plt.show()

List050 = [0, 1, 3, 4, 5, 6, 9, 11, 12, 18, 23, 25, 27, 31, 35, 38, 41, 43, 45, 46, 49, 52, 54, 55, 57, 58, 59, 62, 63, 66, 67, 73, 75, 76, 77, 79, 80, 82, 84, 85, 87, 88, 89, 91, 92, 93, 94, 96, 98, 101, 102, 103, 105, 108, 110, 111, 112, 115, 120, 125, 126, 127, 130, 135, 137, 139, 140, 141, 142, 144, 147, 148, 150, 152, 153, 154, 158, 159, 161, 164, 165, 167, 168, 169, 170, 172, 173, 174, 176, 177, 180, 181, 183, 186, 187, 189, 191, 192, 196, 203, 207, 210, 211, 212, 215, 216, 218, 223, 226, 227, 228, 230, 231, 232, 235, 236, 237, 240, 241, 242, 243, 244, 245, 246, 247, 248, 250, 251, 252, 255, 256, 260, 263, 264, 266, 268, 269, 272, 275, 276, 277, 279, 280, 281, 283, 285, 287, 288, 292, 297, 299, 301, 302, 303, 307, 309, 315, 320, 323, 324, 326, 327, 328, 329, 330, 331, 332, 333, 334, 338, 339, 341, 342, 343, 344, 345, 346, 348, 354, 357, 359, 360, 364, 366, 368, 369, 370, 371, 372, 373, 374, 376, 378, 381, 382, 383, 384, 386, 387, 388, 390, 391, 392, 393, 394, 395, 399, 400, 402, 406, 409, 410, 412, 415, 416, 417, 418, 421, 423, 424, 426, 427, 428, 429, 432, 435, 437, 439, 440, 442, 443, 448, 449, 451, 452, 454, 456, 458, 459, 463, 465, 466, 467, 468, 469, 470, 473, 480, 482, 483, 487, 491, 492, 493, 494, 496, 497, 499]

i_folder = 2

listOfboots = []

for ii in range(500):
    if ii not in List050:
        for i in range(NumOfSeries):
            random.randint(0,LengthOfInterval-1)
    else:
        for i in range(NumOfSeries):
            listOfboots.append(random.randint(i_folder*25,i_folder*25+LengthOfInterval-1))
            
n, bins, patches = plt.hist(listOfboots, num_bins, facecolor = 'green')
plt.show()

List075 = [0, 3, 4, 6, 7, 10, 13, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 27, 29, 32, 33, 35, 37, 40, 41, 42, 43, 46, 48, 49, 52, 55, 57, 58, 59, 61, 62, 64, 67, 70, 71, 72, 76, 77, 78, 79, 80, 83, 84, 86, 88, 89, 93, 94, 96, 97, 98, 99, 100, 104, 105, 108, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 123, 124, 125, 127, 128, 130, 131, 132, 135, 136, 137, 138, 140, 141, 143, 148, 150, 151, 152, 153, 154, 155, 158, 159, 160, 162, 165, 166, 170, 172, 173, 174, 175, 176, 178, 179, 181, 183, 186, 187, 188, 191, 192, 193, 194, 195, 196, 197, 198, 200, 202, 203, 204, 205, 207, 210, 211, 212, 217, 219, 223, 224, 225, 228, 229, 232, 233, 235, 236, 237, 238, 240, 241, 244, 247, 248, 250, 251, 253, 255, 256, 257, 261, 265, 270, 271, 272, 273, 275, 277, 279, 280, 286, 287, 288, 292, 293, 294, 296, 297, 298, 299, 301, 302, 303, 305, 307, 308, 310, 311, 313, 315, 316, 318, 319, 320, 321, 323, 327, 330, 331, 332, 334, 336, 337, 338, 340, 341, 342, 343, 344, 345, 347, 348, 349, 350, 351, 352, 354, 355, 356, 357, 358, 359, 362, 363, 364, 365, 366, 370, 372, 373, 375, 377, 379, 380, 381, 382, 384, 386, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 402, 404, 407, 408, 411, 413, 414, 415, 418, 420, 422, 423, 424, 426, 430, 433, 434, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 449, 451, 453, 454, 456, 458, 459, 460, 461, 462, 463, 465, 466, 467, 468, 469, 470, 471, 472, 473, 475, 476, 477, 482, 483, 485, 486, 489, 490, 491, 493, 496, 498, 499]

i_folder = 3

listOfboots = []

for ii in range(500):
    if ii not in List075:
        for i in range(NumOfSeries):
            random.randint(0,LengthOfInterval-1)
    else:
        for i in range(NumOfSeries):
            listOfboots.append(random.randint(i_folder*25,i_folder*25+LengthOfInterval-1))
            
n, bins, patches = plt.hist(listOfboots, num_bins, facecolor = 'green')
plt.show()

List100 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 19, 20, 22, 24, 26, 27, 31, 32, 33, 34, 36, 38, 39, 42, 44, 45, 46, 47, 48, 49, 51, 52, 53, 55, 56, 57, 58, 60, 61, 62, 63, 64, 67, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 83, 84, 85, 86, 87, 88, 90, 92, 93, 94, 95, 98, 101, 102, 103, 104, 105, 106, 108, 109, 111, 113, 114, 115, 116, 117, 119, 120, 121, 122, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 139, 140, 141, 143, 144, 145, 146, 148, 149, 151, 152, 155, 156, 157, 158, 159, 162, 163, 164, 165, 167, 169, 171, 172, 173, 175, 176, 177, 178, 179, 180, 181, 183, 184, 186, 188, 189, 190, 191, 192, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 205, 206, 207, 208, 210, 211, 214, 215, 216, 217, 220, 221, 222, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 239, 240, 241, 242, 243, 244, 245, 247, 248, 249, 250, 251, 252, 253, 254, 255, 257, 258, 259, 260, 261, 262, 263, 265, 266, 267, 268, 270, 271, 272, 273, 276, 277, 278, 279, 280, 281, 282, 283, 284, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 315, 316, 317, 318, 319, 323, 324, 325, 326, 329, 330, 332, 333, 334, 335, 336, 337, 338, 339, 341, 342, 343, 344, 345, 347, 348, 350, 351, 353, 354, 355, 356, 357, 358, 359, 360, 361, 363, 364, 365, 366, 369, 370, 371, 372, 373, 374, 375, 377, 378, 379, 380, 381, 382, 384, 386, 390, 392, 393, 394, 395, 396, 398, 400, 401, 402, 404, 405, 406, 407, 408, 409, 410, 411, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 425, 427, 428, 429, 431, 432, 433, 434, 435, 436, 438, 439, 442, 444, 445, 446, 449, 450, 451, 452, 454, 455, 457, 458, 459, 460, 461, 463, 465, 466, 468, 469, 470, 471, 472, 473, 475, 476, 477, 479, 480, 481, 482, 483, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497]

i_folder = 4

listOfboots = []

for ii in range(500):
    if ii not in List100:
        for i in range(NumOfSeries):
            random.randint(0,LengthOfInterval-1)
    else:
        for i in range(NumOfSeries):
            listOfboots.append(random.randint(i_folder*25,i_folder*25+LengthOfInterval-1))
            
n, bins, patches = plt.hist(listOfboots, num_bins, facecolor = 'green')
plt.show()

List125 = [0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 51, 52, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 151, 152, 153, 155, 156, 157, 158, 159, 160, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 184, 186, 188, 189, 190, 191, 192, 193, 194, 196, 197, 198, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 215, 216, 218, 219, 220, 221, 222, 223, 224, 226, 227, 228, 229, 230, 231, 232, 233, 234, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 296, 297, 298, 299, 300, 301, 302, 303, 305, 306, 308, 309, 310, 311, 312, 313, 314, 317, 318, 319, 320, 321, 322, 323, 325, 326, 327, 328, 329, 330, 331, 332, 333, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 360, 361, 362, 364, 365, 366, 367, 368, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 393, 394, 395, 396, 397, 398, 399, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 439, 440, 441, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499]

i_folder = 5

listOfboots = []

for ii in range(500):
    if ii not in List125:
        for i in range(NumOfSeries):
            random.randint(0,LengthOfInterval-1)
    else:
        for i in range(NumOfSeries):
            listOfboots.append(random.randint(i_folder*25,i_folder*25+LengthOfInterval-1))
            
n, bins, patches = plt.hist(listOfboots, num_bins, facecolor = 'green')
plt.show()

random.seed(10000)

#6,7,8,9,10,11
List150 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 45, 46, 47, 48, 49, 50, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 69, 70, 71, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 92, 93, 94, 95, 96, 97, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 178, 179, 180, 181, 182, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 199, 201, 202, 203, 204, 205, 206, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 379, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 438, 439, 440, 441, 442, 443, 444, 445, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499]

i_folder = 6

listOfboots = []

for ii in range(500):
    if ii not in List150:
        for i in range(NumOfSeries):
            random.randint(0,LengthOfInterval-1)
    else:
        for i in range(NumOfSeries):
            listOfboots.append(random.randint(i_folder*25,i_folder*25+LengthOfInterval-1))
            
n, bins, patches = plt.hist(listOfboots, num_bins, facecolor = 'green')
plt.show()

List175 = [0, 1, 2, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 51, 52, 53, 54, 55, 57, 59, 60, 61, 62, 63, 64, 65, 67, 69, 70, 71, 72, 73, 75, 76, 77, 78, 79, 80, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 123, 124, 126, 127, 128, 129, 131, 132, 133, 134, 135, 137, 138, 139, 140, 141, 142, 143, 144, 146, 147, 148, 149, 150, 151, 153, 154, 155, 157, 158, 159, 160, 161, 163, 164, 165, 166, 167, 168, 169, 170, 172, 173, 174, 175, 176, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 268, 269, 270, 272, 273, 275, 276, 277, 278, 279, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 311, 312, 313, 315, 316, 317, 318, 319, 322, 323, 324, 325, 326, 327, 328, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 359, 360, 361, 362, 363, 364, 365, 366, 367, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 381, 382, 383, 384, 385, 386, 388, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415, 417, 418, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 457, 459, 460, 461, 462, 464, 465, 466, 467, 469, 470, 471, 473, 474, 475, 476, 477, 479, 480, 481, 482, 483, 485, 487, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498]

i_folder = 7

listOfboots = []

for ii in range(500):
    if ii not in List175:
        for i in range(NumOfSeries):
            random.randint(0,LengthOfInterval-1)
    else:
        for i in range(NumOfSeries):
            listOfboots.append(random.randint(i_folder*25,i_folder*25+LengthOfInterval-1))
            
n, bins, patches = plt.hist(listOfboots, num_bins, facecolor = 'green')
plt.show()

List200 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 15, 16, 17, 19, 20, 21, 23, 24, 25, 27, 28, 29, 30, 31, 33, 36, 37, 38, 39, 40, 43, 45, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 87, 88, 89, 90, 91, 93, 94, 95, 96, 97, 98, 99, 100, 101, 103, 104, 105, 106, 107, 109, 111, 112, 113, 114, 115, 116, 117, 119, 120, 121, 122, 123, 125, 126, 127, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 190, 191, 192, 193, 194, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 208, 209, 210, 211, 212, 213, 214, 215, 216, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 235, 237, 239, 240, 241, 242, 243, 244, 245, 246, 247, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 266, 267, 268, 269, 270, 271, 272, 274, 275, 276, 277, 278, 279, 280, 281, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303, 304, 306, 307, 308, 309, 310, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 334, 335, 336, 337, 338, 339, 340, 341, 342, 345, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 362, 363, 364, 366, 367, 368, 370, 371, 372, 373, 374, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 387, 388, 389, 390, 391, 392, 394, 395, 396, 397, 398, 400, 401, 402, 403, 404, 405, 407, 408, 409, 411, 412, 413, 414, 415, 416, 418, 419, 420, 421, 422, 424, 425, 426, 427, 428, 429, 430, 432, 433, 434, 435, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 455, 456, 457, 458, 460, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 474, 475, 476, 477, 478, 479, 481, 483, 484, 485, 486, 487, 489, 490, 492, 495, 496, 497, 498, 499]

i_folder = 8

listOfboots = []

for ii in range(500):
    if ii not in List200:
        for i in range(NumOfSeries):
            random.randint(0,LengthOfInterval-1)
    else:
        for i in range(NumOfSeries):
            listOfboots.append(random.randint(i_folder*25,i_folder*25+LengthOfInterval-1))
            
n, bins, patches = plt.hist(listOfboots, num_bins, facecolor = 'green')
plt.show()

List225 = [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 16, 18, 19, 20, 21, 22, 24, 25, 28, 29, 30, 31, 32, 33, 34, 35, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 71, 73, 74, 75, 76, 77, 78, 79, 80, 81, 83, 84, 86, 87, 88, 89, 90, 92, 93, 95, 97, 98, 99, 100, 101, 102, 103, 104, 105, 107, 108, 110, 111, 112, 113, 114, 115, 116, 117, 118, 122, 123, 124, 125, 126, 127, 128, 129, 130, 134, 135, 136, 138, 139, 140, 141, 143, 144, 146, 148, 149, 150, 151, 152, 153, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 194, 195, 196, 197, 198, 199, 200, 202, 203, 204, 206, 207, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 253, 254, 256, 257, 258, 259, 262, 263, 264, 266, 267, 268, 269, 271, 272, 273, 274, 276, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 291, 292, 293, 295, 296, 297, 298, 299, 300, 301, 302, 304, 305, 306, 307, 308, 309, 310, 312, 313, 314, 315, 316, 318, 319, 320, 321, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 347, 348, 349, 350, 351, 353, 354, 355, 356, 357, 359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 386, 387, 388, 389, 390, 391, 392, 394, 395, 396, 397, 398, 399, 400, 401, 402, 405, 406, 407, 408, 410, 411, 412, 414, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 442, 443, 444, 445, 446, 447, 448, 450, 451, 452, 453, 454, 455, 456, 457, 458, 459, 460, 461, 463, 464, 467, 468, 469, 470, 471, 472, 474, 475, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498]

i_folder = 9

listOfboots = []

for ii in range(500):
    if ii not in List225:
        for i in range(NumOfSeries):
            random.randint(0,LengthOfInterval-1)
    else:
        for i in range(NumOfSeries):
            listOfboots.append(random.randint(i_folder*25,i_folder*25+LengthOfInterval-1))
            
n, bins, patches = plt.hist(listOfboots, num_bins, facecolor = 'green')
plt.show()

List250 = [0, 1, 4, 6, 7, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 38, 39, 41, 42, 43, 46, 47, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 60, 61, 62, 63, 64, 65, 66, 68, 69, 71, 74, 75, 77, 78, 80, 81, 82, 83, 84, 85, 87, 88, 89, 92, 93, 95, 96, 98, 99, 100, 102, 103, 104, 105, 107, 108, 109, 110, 112, 114, 115, 116, 118, 119, 120, 121, 122, 125, 127, 129, 131, 132, 133, 134, 135, 136, 137, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 151, 152, 153, 154, 156, 157, 158, 159, 161, 163, 164, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 184, 186, 188, 189, 190, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 211, 213, 214, 216, 217, 219, 220, 221, 222, 223, 224, 225, 226, 227, 229, 230, 231, 232, 233, 234, 235, 236, 237, 239, 240, 241, 242, 243, 244, 245, 246, 248, 249, 250, 251, 252, 254, 255, 256, 257, 258, 259, 260, 262, 263, 265, 266, 268, 270, 272, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 285, 286, 287, 289, 290, 292, 295, 296, 298, 299, 300, 301, 302, 303, 304, 306, 307, 308, 310, 311, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 327, 328, 330, 331, 333, 335, 339, 340, 342, 343, 344, 347, 348, 349, 350, 351, 352, 353, 357, 359, 361, 362, 363, 364, 365, 367, 368, 369, 370, 371, 372, 373, 374, 375, 377, 378, 380, 381, 383, 385, 386, 387, 391, 392, 394, 395, 396, 398, 399, 400, 401, 402, 404, 405, 406, 407, 408, 409, 410, 411, 413, 414, 415, 417, 419, 421, 422, 423, 424, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435, 437, 438, 439, 440, 441, 442, 443, 445, 447, 448, 449, 450, 451, 452, 453, 454, 455, 456, 459, 460, 462, 463, 464, 466, 467, 468, 469, 470, 471, 472, 474, 475, 476, 477, 478, 479, 480, 481, 484, 485, 486, 488, 489, 490, 492, 493, 495, 496, 497, 498, 499]

i_folder = 10

listOfboots = []

for ii in range(500):
    if ii not in List250:
        for i in range(NumOfSeries):
            random.randint(0,LengthOfInterval-1)
    else:
        for i in range(NumOfSeries):
            listOfboots.append(random.randint(i_folder*25,i_folder*25+LengthOfInterval-1))
            
n, bins, patches = plt.hist(listOfboots, num_bins, facecolor = 'green')
plt.show()

List275 = [2, 4, 5, 7, 13, 15, 19, 20, 22, 25, 26, 29, 31, 32, 33, 34, 41, 42, 43, 44, 46, 47, 48, 49, 52, 54, 55, 56, 57, 59, 60, 61, 62, 63, 65, 67, 70, 72, 73, 78, 79, 80, 81, 82, 83, 84, 87, 88, 89, 91, 93, 94, 95, 96, 97, 98, 99, 102, 108, 109, 112, 113, 117, 118, 123, 124, 127, 129, 131, 133, 135, 136, 138, 141, 145, 147, 152, 154, 156, 160, 162, 163, 165, 171, 173, 174, 175, 179, 182, 183, 185, 186, 188, 189, 190, 196, 197, 206, 207, 208, 210, 211, 212, 213, 215, 216, 218, 219, 220, 221, 222, 225, 226, 227, 237, 238, 239, 240, 246, 247, 249, 250, 260, 262, 263, 265, 266, 267, 268, 271, 275, 276, 281, 283, 284, 288, 289, 290, 291, 292, 293, 294, 295, 297, 298, 303, 305, 306, 309, 310, 313, 314, 316, 317, 318, 321, 322, 323, 325, 326, 327, 328, 329, 331, 332, 334, 335, 338, 342, 343, 344, 345, 346, 349, 350, 351, 353, 356, 357, 360, 363, 366, 367, 369, 370, 371, 373, 374, 377, 380, 381, 384, 385, 388, 389, 391, 393, 397, 398, 399, 400, 402, 403, 405, 407, 409, 410, 413, 414, 417, 419, 420, 424, 426, 430, 431, 432, 435, 437, 441, 442, 443, 444, 446, 447, 448, 450, 451, 452, 453, 454, 455, 457, 461, 462, 465, 470, 471, 472, 480, 483, 484, 485, 487, 488, 489, 490, 492, 494, 496, 497, 498, 499]

i_folder = 11

listOfboots = []

for ii in range(500):
    if ii not in List275:
        for i in range(NumOfSeries):
            random.randint(0,LengthOfInterval-1)
    else:
        for i in range(NumOfSeries):
            listOfboots.append(random.randint(i_folder*25,i_folder*25+LengthOfInterval-1))
            
n, bins, patches = plt.hist(listOfboots, num_bins, facecolor = 'green')
plt.show()

random.seed(10000)

#12,13,14,15
List300 = [1, 2, 4, 6, 8, 13, 14, 17, 19, 26, 27, 28, 29, 30, 33, 36, 39, 42, 45, 47, 48, 54, 59, 60, 61, 62, 64, 69, 70, 71, 76, 81, 87, 91, 92, 94, 96, 100, 103, 105, 107, 108, 113, 114, 116, 117, 118, 122, 123, 124, 125, 128, 130, 131, 133, 134, 138, 147, 151, 155, 156, 168, 170, 171, 172, 177, 183, 184, 187, 188, 189, 190, 191, 194, 195, 201, 203, 205, 208, 209, 215, 216, 218, 219, 223, 226, 229, 230, 231, 233, 238, 240, 243, 245, 246, 254, 255, 264, 268, 270, 272, 276, 278, 282, 293, 298, 309, 311, 318, 319, 323, 326, 329, 340, 345, 346, 347, 349, 351, 352, 353, 358, 366, 371, 372, 381, 386, 389, 394, 399, 402, 403, 406, 408, 413, 414, 415, 416, 422, 426, 428, 431, 438, 439, 440, 443, 446, 453, 462, 470, 472, 474, 477, 483, 484, 486, 492, 494, 495]

i_folder = 12

listOfboots = []

for ii in range(500):
    if ii not in List300:
        for i in range(NumOfSeries):
            random.randint(0,LengthOfInterval-1)
    else:
        for i in range(NumOfSeries):
            listOfboots.append(random.randint(i_folder*25,i_folder*25+LengthOfInterval-1))
            
n, bins, patches = plt.hist(listOfboots, num_bins, facecolor = 'green')
plt.show()

List325 = [2, 4, 5, 6, 7, 12, 18, 19, 20, 23, 29, 30, 31, 33, 40, 45, 48, 54, 56, 57, 59, 64, 65, 68, 71, 72, 74, 77, 80, 84, 86, 87, 89, 90, 97, 103, 107, 108, 109, 110, 112, 114, 115, 116, 118, 119, 120, 123, 127, 128, 129, 131, 138, 139, 140, 141, 142, 145, 152, 153, 158, 165, 172, 175, 176, 177, 178, 180, 181, 185, 187, 188, 197, 199, 204, 205, 208, 209, 210, 217, 219, 221, 224, 226, 227, 229, 231, 235, 236, 237, 238, 244, 245, 255, 256, 257, 259, 260, 262, 263, 265, 266, 267, 271, 275, 279, 283, 284, 286, 288, 290, 296, 299, 301, 303, 304, 305, 308, 310, 311, 312, 313, 315, 319, 320, 325, 330, 333, 338, 349, 354, 355, 356, 360, 362, 364, 366, 368, 375, 376, 377, 379, 381, 383, 385, 387, 392, 395, 399, 401, 405, 407, 410, 411, 415, 417, 418, 419, 426, 427, 428, 436, 437, 439, 450, 452, 454, 460, 466, 468, 473, 474, 476, 477, 479, 480, 484, 485, 488, 490, 491, 495, 496, 498, 499]

i_folder = 13

listOfboots = []

for ii in range(500):
    if ii not in List325:
        for i in range(NumOfSeries):
            random.randint(0,LengthOfInterval-1)
    else:
        for i in range(NumOfSeries):
            listOfboots.append(random.randint(i_folder*25,i_folder*25+LengthOfInterval-1))
            
n, bins, patches = plt.hist(listOfboots, num_bins, facecolor = 'green')
plt.show()

List350 = [4, 8, 9, 11, 12, 15, 18, 21, 22, 24, 30, 35, 39, 63, 65, 68, 70, 76, 77, 78, 83, 86, 101, 102, 103, 110, 113, 120, 122, 125, 127, 135, 139, 144, 145, 163, 166, 168, 183, 187, 207, 212, 219, 226, 232, 239, 243, 248, 255, 258, 262, 271, 272, 274, 276, 283, 293, 303, 305, 326, 329, 339, 346, 352, 356, 364, 370, 375, 377, 382, 386, 389, 394, 404, 405, 407, 410, 412, 413, 423, 433, 436, 451, 454, 456, 457, 459, 460, 464, 465, 477, 487, 489, 494, 495, 499]

i_folder = 14

listOfboots = []

for ii in range(500):
    if ii not in List350:
        for i in range(NumOfSeries):
            random.randint(0,LengthOfInterval-1)
    else:
        for i in range(NumOfSeries):
            listOfboots.append(random.randint(i_folder*25,i_folder*25+LengthOfInterval-1))
            
n, bins, patches = plt.hist(listOfboots, num_bins, facecolor = 'green')
plt.show()

List375 = [2, 3, 13, 14, 24, 45, 73, 77, 92, 120, 125, 127, 130, 136, 163, 164, 172, 175, 188, 190, 207, 212, 216, 219, 231, 244, 262, 267, 269, 278, 285, 290, 301, 307, 309, 315, 317, 321, 325, 326, 328, 331, 342, 355, 360, 363, 375, 378, 383, 387, 396, 425, 444, 446, 462, 474, 476, 477, 496, 497]

i_folder = 15

listOfboots = []

for ii in range(500):
    if ii not in List375:
        for i in range(NumOfSeries):
            random.randint(0,LengthOfInterval-1)
    else:
        for i in range(NumOfSeries):
            listOfboots.append(random.randint(i_folder*25,i_folder*25+LengthOfInterval-1))
            
n, bins, patches = plt.hist(listOfboots, num_bins, facecolor = 'green')
plt.show()

random.seed(10000)

#-1
List400 = [12, 24, 31, 39, 66, 86, 87, 93, 108, 109, 120, 133, 136, 151, 153, 171, 183, 186, 188, 190, 193, 203, 204, 207, 222, 227, 231, 264, 297, 306, 307, 311, 328, 333, 336, 338, 351, 353, 356, 360, 361, 362, 364, 365, 389, 393, 414, 416, 427, 428, 436, 445, 446, 455, 456, 470, 483, 488, 493, 495, 496]

i_folder = 16

listOfboots = []

for ii in range(500):
    if ii not in List400:
        for i in range(NumOfSeries):
            random.randint(0,LengthOfInterval-1)
    else:
        for i in range(NumOfSeries):
            listOfboots.append(random.randint(i_folder*25,i_folder*25+LengthOfInterval-1))
            
n, bins, patches = plt.hist(listOfboots, num_bins, facecolor = 'green')
plt.show()