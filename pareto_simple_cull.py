# downloaded from http://pythonfiddle.com/pareto-simple-cull/
# this function finds which points lies on pareto frontier
# input format: inputPoints: list of points, where items are distances from thresholds,
#                            eg [[dist1, dist2, ...], [dist1, dist2, ...]]
#               dominates: function which decides which point is dominated
# output format: return indexes in list of non-dominated points
def simple_cull(inputPoints, dominates):

    paretoPoints = list()
    candidateRowNr = 0
    index = 0

    while len(inputPoints):
        candidateRow = inputPoints[candidateRowNr]
        inputPoints.remove(candidateRow)
        rowNr = 0
        nonDominated = True

        while len(inputPoints) != 0 and rowNr < len(inputPoints):
            row = inputPoints[rowNr]
            if dominates(candidateRow, row):
                inputPoints.remove(row)
            elif dominates(row, candidateRow):
                nonDominated = False
                rowNr += 1
            else:
                rowNr += 1

        if nonDominated:
            paretoPoints.append(index+1)
        index += 1

    return paretoPoints


def dominates_min(row, anotherRow):
    return sum([row[x] <= anotherRow[x] for x in range(len(row))]) == len(row)  # minimization domination


# return simple distance from threshold
# this function is used in process_predictions.py, numpy version use its own vectorized version of this function
# it produces value which represents distance from threshold
# input format: point: float number, predicted value
#               threshold: list with threshold, e.g. ['more', -3.2] or ['between', -0.5, 2.8]
# output format: return 0 if the value matches the threshold,
#                if not it returns value which represents distance from threshold
def get_distance_from_threshold(point, threshold):
    if threshold[0] == 'more':
        return 0 if point > threshold[1] else threshold[1] - point
    elif threshold[0] == 'less':
        return 0 if point < threshold[1] else point - threshold[1]
    else:
        if point >= threshold[1] and point <= threshold[2]:
            return 0
        else:
            return point - threshold[2] if point > threshold[2] else threshold[1] - point

