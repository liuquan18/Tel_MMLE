# print hellow world
print("Hello World")

# count the extreme events from a list
# the extremes are defined as the values that are greater than 0.5 (positive)
# or less than -0.5 (negative)
def count_extremes(list, gt):
    count = 0
    for i in range(len(list)):
        if gt:
            if list[i] > 0.5:
                count += 1
        else:
            if list[i] < -0.5:
                count += 1
    return count

# count the number of positive and negative extremes
list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
ext_pos = count_extremes(list, True)
ext_neg = count_extremes(list, False)