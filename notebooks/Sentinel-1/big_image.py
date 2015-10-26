__author__ = 'Denis Spiridonov'

import os

from pylab import *
from PIL import Image
from math import ceil

def save_image_part(array, x, y, out_dir, vmin, vmax):
    filepath = os.path.join(out_dir, '%d_%d.png' % (x, y))
    imsave(filepath, array, vmin=vmin, vmax=vmax)

def save_big_image(source_array, result_image_path, vmin, vmax,
                   part_size=4000):
    shape = source_array.shape
    temp_dir = '/tmp/save_image_sentinel'

    if not os.path.isdir(temp_dir):
        os.mkdir(temp_dir)

    print "Start create images"
    part_x_count = int(ceil(source_array.shape[0]/float(part_size)))
    part_y_count = int(ceil(source_array.shape[1]/float(part_size)))

    for i in range(part_x_count):
        temp_x = source_array[part_size*i:part_size*(i+1), :]
        for j in range(part_y_count):
#             print "%d:%d" % (part_size*j, part_size*(j+1))
#             temp_y = temp_x[:, part_size*j:part_size*(j+1)]
            save_image_part(temp_x[:, part_size*j:part_size*(j+1)], i, j,
                            temp_dir, vmin, vmax)

    print "Start paste to big image"
    if os.path.isfile(result_image_path):
        os.remove(result_image_path)

    img = Image.new('RGBA', (shape[1], shape[0]))

    for i in range(part_x_count):
        for j in range(part_y_count):
            temp_file_path = os.path.join(temp_dir, '%d_%d.png' % (i, j))
            img_part = Image.open(temp_file_path)
            img.paste(img_part, (part_size*j, part_size*i))
            os.remove(temp_file_path)

    png_info = img_part.info
    img.save(result_image_path, **png_info)
    print "Image saved successfully"

# save_big_image(result, '/tmp/save_image', 0, 2, 4000)
