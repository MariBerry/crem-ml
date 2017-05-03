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

