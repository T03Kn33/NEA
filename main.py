import time as time
#import tkinter as tk
#window = tk.Tk()
#window.title("Hello world")
#window.geometry("300x300")
#
#hello = tk.Label(text="Hello world!")
#hello.pack()
#button = tk.Button(text="Click me!")
#button.pack()
#tk.mainloop()


def PhysicsExplanation():
    print(
        "\nThe permittivity of a material is the ratio of the permeability of the"
        " material to the permeability of vacuum.(This is a placeholder)")


#GeeksforGeeks | Source:(https://www.geeksforgeeks.org/python-program-for-merge-sort/)
def merge(arr, l, m, r):
    n1 = m - l + 1
    n2 = r - m

    # create temp arrays
    L = [0] * (n1)
    R = [0] * (n2)

    # Copy data to temp arrays L[] and R[]
    for i in range(0, n1):
        L[i] = arr[l + i]

    for j in range(0, n2):
        R[j] = arr[m + 1 + j]

    # Merge the temp arrays back into arr[l..r]
    i = 0  # Initial index of first subarray
    j = 0  # Initial index of second subarray
    k = l  # Initial index of merged subarray

    while i < n1 and j < n2:
        if L[i] <= R[j]:
            arr[k] = L[i]
            i += 1
        else:
            arr[k] = R[j]
            j += 1
        k += 1

    # Copy the remaining elements of L[], if there
    # are any
    while i < n1:
        arr[k] = L[i]
        i += 1
        k += 1

    # Copy the remaining elements of R[], if there
    # are any
    while j < n2:
        arr[k] = R[j]
        j += 1
        k += 1


# l is for left index and r is right index of the
# sub-array of arr to be sorted


def mergeSort(arr, l, r):
    if l < r:

        # Same as (l+r)//2, but avoids overflow for
        # large l and h
        m = l + (r - l) // 2

        # Sort first and second halves
        mergeSort(arr, l, m)
        mergeSort(arr, m + 1, r)
        merge(arr, l, m, r)


#GeeksforGeeks | Source:(https://www.geeksforgeeks.org/python-program-for-binary-search/)
# Returns index of x in arr if present, else -1
def binary_search(arr, low, high, x):

    # Check base case
    if high >= low:

        mid = (high + low) // 2

        # If element is present at the middle itself
        if arr[mid] == x:
            return mid

        # If element is smaller than mid, then it can only
        # be present in left subarray
        elif arr[mid] > x:
            return binary_search(arr, low, mid - 1, x)

        # Else the element can only be present in right subarray
        else:
            return binary_search(arr, mid + 1, high, x)

    else:
        # Element is not present in the array
        return -1


#Unsorted list of permittivities from https://www.engineeringtoolbox.com/relative-permittivity-d_1660.html
permittivity = [
    4.5, 5, 3.6, 3, 1.5, 83.6, 6.2, 2.3, 37, 2.2, 33.6, 15.8, 3.5, 4.7, 21.8,
    18.9, 8, 3.9, 3.2, 22, 6, 24.3, 16.5, 84, 2.5, 4, 42, 4.8, 7.8, 7.3, 60,
    6.7, 52, 38, 21, 3.1, 2.6, 4.2, 1.4, 28, 17.7, 20.7, 2.4, 25, 42.5, 23.1,
    4.1, 5.6, 2.1, 3.4, 2.8, 5.5, 1, 2
]

relativePermittivities = {
    "product0": {
        "name": "Vacuum",
        "relative_permittivity": 1
    },
    "product1": {
        "name": "Air, Liquid (-191 oC)",
        "relative_permittivity": 1.4
    },
    "product2": {
        "name": "Insulation of telephone cables",
        "relative_permittivity": 1.5
    },
    "product3": {
        "name": "Polytetrafluoroethylene (PTFE)",
        "relative_permittivity": 2
    },
    "product4": {
        "name": "Teflon",
        "relative_permittivity": 2.1
    },
    "product5": {
        "name": "Petroleum",
        "relative_permittivity": 2.2
    },
    "product6": {
        "name": "Paper",
        "relative_permittivity": 2.3
    },
    "product7": {
        "name": "R12 Dichlorodifluoromethane (70°F)",
        "relative_permittivity": 2.4
    },
    "product8": {
        "name": "Casting compound",
        "relative_permittivity": 2.5
    },
    "product9": {
        "name": "Carbon disulfide",
        "relative_permittivity": 2.6
    },
    "product10": {
        "name": "Polyamide",
        "relative_permittivity": 2.8
    },
    "product11": {
        "name": "Polystyrene",
        "relative_permittivity": 3
    },
    "product12": {
        "name": "Mylar",
        "relative_permittivity": 3.1
    },
    "product13": {
        "name": "Plexiglass",
        "relative_permittivity": 3.2
    },
    "product14": {
        "name": "Polymide",
        "relative_permittivity": 3.4
    },
    "product15": {
        "name": "Sulfur",
        "relative_permittivity": 3.5
    },
    "product16": {
        "name": "Araldite",
        "relative_permittivity": 3.6
    },
    "product17": {
        "name": "Silicon dioxide",
        "relative_permittivity": 3.9
    },
    "product18": {
        "name": "Slate",
        "relative_permittivity": 4
    },
    "product19": {
        "name": "Acetic Acid (36°F)",
        "relative_permittivity": 4.1
    },
    "product20": {
        "name": "Insulation of high voltage cables",
        "relative_permittivity": 4.2
    },
    "product21": {
        "name": "Concrete, quartz",
        "relative_permittivity": 4.5
    },
    "product22": {
        "name": "Caster oil",
        "relative_permittivity": 4.7
    },
    "product23": {
        "name": "Chloroform (68°F)",
        "relative_permittivity": 4.8
    },
    "product24": {
        "name": "Paper, impregnated",
        "relative_permittivity": 5
    },
    "product25": {
        "name": "Aniline (212°F)",
        "relative_permittivity": 5.5
    },
    "product26": {
        "name": "Sodium chloride",
        "relative_permittivity": 5.6
    },
    "product27": {
        "name": "Ethyl Acetate (77°F)",
        "relative_permittivity": 6
    },
    "product28": {
        "name": "Acetic Acid (68°F)",
        "relative_permittivity": 6.2
    },
    "product29": {
        "name": "Neoprene",
        "relative_permittivity": 6.7
    },
    "product30": {
        "name": "Aniline (68°F)",
        "relative_permittivity": 7.3
    },
    "product31": {
        "name": "Aniline (32°F)",
        "relative_permittivity": 7.8
    },
    "product32": {
        "name": "Body tissue, Marble",
        "relative_permittivity": 8
    },
    "product33": {
        "name": "Acetyl Chloride (68°F)",
        "relative_permittivity": 15.8
    },
    "product34": {
        "name": "Ammonia (69°F)",
        "relative_permittivity": 16.5
    },
    "product35": {
        "name": "Acetone (127°F)",
        "relative_permittivity": 17.7
    },
    "product36": {
        "name": "Ammonia (40°F)",
        "relative_permittivity": 18.9
    },
    "product37": {
        "name": "Acetone (77°F)",
        "relative_permittivity": 20.7
    },
    "product38": {
        "name": "Ceramic, MgNb2O6",
        "relative_permittivity": 21
    },
    "product39": {
        "name": "Acetaldehyde (41°F)",
        "relative_permittivity": 21.8
    },
    "product40": {
        "name": "Ammonia (-30°F)",
        "relative_permittivity": 22
    },
    "product41": {
        "name": "Acetyl Acetone (68°F)",
        "relative_permittivity": 23.1
    },
    "product42": {
        "name": "Ethanol (25 oC, 77°F)",
        "relative_permittivity": 24.3
    },
    "product43": {
        "name": "Ammonia (-74°F), Ceramic, ZnNb2O6",
        "relative_permittivity": 25
    },
    "product44": {
        "name": "Ceramic, MgTa2O6",
        "relative_permittivity": 28
    },
    "product45": {
        "name": "Methanol (20 oC))",
        "relative_permittivity": 33.6
    },
    "product46": {
        "name": "Ethylene glycol",
        "relative_permittivity": 37
    },
    "product47": {
        "name": "Ceramic, ZnTa2O6",
        "relative_permittivity": 38
    },
    "product48": {
        "name": "Furfural",
        "relative_permittivity": 42
    },
    "product49": {
        "name": "Glycerol (77°F)",
        "relative_permittivity": 42.5
    },
    "product50": {
        "name": "Hydrazine (20 oC)",
        "relative_permittivity": 52
    },
    "product51": {
        "name": "Hydrogen peroxide (25 oC)",
        "relative_permittivity": 60
    },
    "product52": {
        "name": "Hydrofluoric acid (0 oC)",
        "relative_permittivity": 83.6
    },
    "product53": {
        "name": "Formamide (20 oC), Sulfuric acid (20 oC)",
        "relative_permittivity": 84
    },
}  #Dictionary containing all the permittivities and associated materials

n = len(permittivity)

mergeSort(permittivity, 0, n - 1)
print("\nSorted array is")
for i in range(n - 1):
    print(permittivity[i], end=", ")
print(permittivity[-1])


def permittivityValidation():
    try:
        x = float(input("\nEnter a permittivity to search for: "))
        if type(x) is int or type(x) is float:
            return x

        else:
            print("Please enter a number")
            return permittivityValidation()
    except TypeError:
        print("Please enter a number")
        return permittivityValidation()

    except ValueError:
        print("Please enter a number")
        return permittivityValidation()


x = permittivityValidation()


# Function call
def findPermittivity(permittivity, x):
    result = binary_search(permittivity, 0, len(permittivity) - 1, x)

    if result != -1:
        print("Element is present at index", str(result))
        return result
    else:
        print("Element is not present in array")
        return -1


permittivityIndex = findPermittivity(permittivity, x)


def commonPermittivity():
    if permittivityIndex != -1:
        print("\nThe permittivity of the dielectric is the same as " +
              relativePermittivities[f"product{permittivityIndex}"]["name"])


def findCapacitance(ε, area, distance):
    return (ε * area) / distance


commonPermittivity()

u = 0.0  #U = Initial velocity
pd = 10.0  #pd = Potential difference
gap = 0.5  #gap = Gap between the two plates
ε0 = 8.854187817 * (10**(-12))  #ε0 = Permittivity of free space
εr = 1.0  #εr = Relative permittivity of the dielectric
mass = 0.05  #mass = Mass of the ball
charge = 1.0 * (10**(-6))  #charge = Charge of the ball
rho = 100.0  #rho = Density of the dielectric
cd = 0.47  #cd = Drag coefficient
pi = 3.14159
r = 0.05  #r = Radius of the ball
t = 0.1  #t = Time interval
areaPlates = 0.01
drag = 0


class baselineValues:  #Class for the baseline values of the componenents

    def __init__(self, u, pd, gap, ε0, εr, mass, charge, rho, cd, pi, r, t,
                 areaPlates, drag):
        self.u = u
        self.pd = pd
        self.gap = gap
        self.ε0 = ε0
        self.εr = εr
        self.mass = mass
        self.charge = charge
        self.rho = rho
        self.cd = cd
        self.pi = pi
        self.r = r
        self.crossArea = r * r * pi
        self.t = t
        self.areaPlates = areaPlates
        self.drag = drag

        def getU(self):
            return self.u

        def getPD(self):
            return self.pd

        def getGap(self):
            return self.gap

        def getε0(self):
            return self.ε0

        def getεr(self):
            return self.εr

        def getMass(self):
            return self.mass

        def getCharge(self):
            return self.charge

        def getRho(self):
            return self.rho

        def getCD(self):
            return self.cd

        def getPi(self):
            return self.pi

        def getR(self):
            return self.r

        def getCrossArea(self):
            return self.crossArea

        def getT(self):
            return self.t

        def getDrag(self):
            return self.drag


#When drag is close to force, make them equal to each other
def dragCheck(drag, force):
    if abs((drag - force) / force) >= 0.05:
        drag = force


class Ball:
    global charge
    global mass

    def __init__(self, u, v, mass, acceleration, charge, r, t):
        self.u = u
        self.v = v
        self.mass = mass
        self.acceleration = acceleration
        self.t = t
        self.s = (self.u * self.t) + (0.5 * self.acceleration * self.t *
                                      self.t)
        self.charge = charge
        self.r = r
        self.crossArea = r * r * pi
        self.drag = -(0.5 * dielectric.cd * dielectric.rho * self.crossArea *
                      self.v * self.v)  #Equation for drag

    def getU(self):
        return self.u

    def getV(self):
        return self.v

    def getMass(self):
        return self.mass

    def getAcceleration(self):
        return self.acceleration

    def getS(self):
        return self.s

    def getCharge(self):
        return self.charge

    def getDrag(self):
        return self.drag

    def getR(self):
        return self.r

    def getcrossArea(self):
        return self.crossArea

    def updateU(self, u):
        self.u = u

    def updateV(self):
        self.v = self.u + self.acceleration * self.t

    def updateMass(self, mass):
        self.mass = mass

    def updateAcceleration(self):
        self.acceleration = (field.force + self.drag) / self.mass

    def updateS(self):
        self.s = self.u * self.t + (0.5 * self.acceleration * self.t * self.t)

    def updateCharge(self, charge):
        self.charge = charge

    def updateDrag(self, drag):
        self.drag = drag

    def updateR(self, r):
        self.r = r

    def reverseU(self):
        self.u = -self.u

    def reverseV(self):
        self.v = -self.v

    def reverseAcceleration(self):
        self.acceleration = -self.acceleration

    def reverseDrag(self):
        self.drag = -self.drag

    def classVariables(self):
        print("u = " + str(self.u))
        print("v = " + str(self.v))
        print("mass = " + str(self.mass))
        print("acceleration = " + str(self.acceleration))
        print("s = " + str(self.s))
        print("charge = " + str(self.charge))
        print("r = " + str(self.r))
        print("crossArea = " + str(self.crossArea))
        print("drag = " + str(self.drag))


class Plates:
    global pd
    global gap
    global areaPlates

    def __init__(self, pd, gap, areaPlates):
        self.pd = pd
        self.gap = gap
        self.areaPlates = areaPlates

    def getPD(self):
        return self.pd

    def getGapPlates(self):
        return self.gap

    def getAreaPlates(self):
        return self.areaPlates

    def updatePD(self, pd):
        self.pd = pd

    def updateGap(self, gap):
        self.gap = gap

    def updateAreaPlates(self, areaPlates):
        self.areaPlates = areaPlates

    def classVariables(self):
        print("pd = " + str(self.pd))
        print("gap = " + str(self.gap))
        print("areaPlates = " + str(self.areaPlates))


#Try use inheritance for field and dielectric | (https://www.w3schools.com/python/python_inheritance.asp)
class Dielectric(Plates):
    global ε
    global rho

    def __init__(self, pd, gap, areaPlates, ε0, εr, rho, cd):
        super().__init__(pd, gap, areaPlates)
        self.ε0 = ε0
        self.εr = εr
        self.ε = self.ε0 * self.εr
        self.rho = rho
        self.cd = cd

    def getε0(self):
        return self.ε0

    def getεr(self):
        return self.εr

    def getε(self):
        return self.ε

    def getRho(self):
        return self.rho

    def getCD(self):
        return self.cd

    def updateε0(self, ε0):
        self.ε0 = ε0

    def updateεr(self, εr):
        self.εr = εr

    def updateε(self, ε):
        self.ε = ε

    def updateRho(self, rho):
        self.rho = rho

    def updateCD(self, cd):
        self.cd = cd

    def classVariables(self):
        print("ε0 = " + str(self.ε0))
        print("εr = " + str(self.εr))
        print("ε = " + str(self.ε))
        print("rho = " + str(self.rho))


class Field(Plates):
    global force
    global fieldStrength

    def __init__(self, pd, gap, areaPlates, fieldStrength, force):
        super().__init__(pd, gap, areaPlates)
        self.fieldStrength = fieldStrength
        self.force = force

    def getFieldStrength(self):
        return self.fieldStrength

    def getForce(self):
        return self.force

    def updateFieldStrength(self, fieldStrength):
        self.fieldStrength = fieldStrength

    def updateForce(self, force):
        self.force = force

    def reverseForce(self):
        self.force = -self.force

    def classVariables(self):
        print("fieldStrength = " + str(self.fieldStrength))
        print("force = " + str(self.force))


baseline = baselineValues(u, pd, gap, ε0, εr, mass, charge, rho, cd, pi, r, t,
                          areaPlates, drag)

plates = Plates(baseline.pd, baseline.gap, baseline.areaPlates)
dielectric = Dielectric(plates.pd, plates.gap, plates.areaPlates, baseline.ε0,
                        baseline.εr, baseline.rho, baseline.cd)
field = Field(plates.pd, plates.gap, plates.areaPlates,
              (plates.pd / plates.gap),
              baseline.charge * (plates.pd / plates.gap))
ball = Ball(baseline.u, 0, baseline.mass,
            (field.force + baseline.drag) / baseline.mass, baseline.charge,
            baseline.r, baseline.t)

#(plates.pd/plates.gap) = fieldStrength | ball.charge * (plates.pd / plates.gap) = force
for i in range(50):
    print(f"Iteration {i}")
    print("\n")
    print(ball.classVariables())
    print(f"field strength = {field.force}")
    ball.updateAcceleration()
    ball.updateV()
    ball.updateS()
    ball.updateU(ball.v)
    ball.updateDrag(dielectric.cd * dielectric.rho * ball.crossArea *
                    ball.v * ball.v)
    dragCheck(ball.drag, field.force)
    
