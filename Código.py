import matplotlib.pyplot as plt
import numpy as np

def dataFH(path):
    
    U1, Ia, minimos = [], [], []
    
    file = open(path, "r", encoding="utf8")
    
    file.readline(), file.readline(), file.readline()
    
    for linea in file:
        
        l = linea.split()
        l0 = float(l[0].replace(",", "."))
        l1 = float(l[1].replace(",", "."))
        
        U1.append(l0), Ia.append(l1)
        
        try:
            
            l3 = l[2]
            minimos.append(l0)
            
        except:
            
            continue
        
    file.close()
    return U1, Ia, minimos

def plotIa(path, save=False, figs=(9, 6), title="Gráfico de 3V"):
    
    U1, Ia, minimos = dataFH(path)
    
    plt.figure(figsize=(figs))
    plt.plot(U1, Ia)
    plt.xlabel("$U_1 (V)$")
    plt.ylabel("$I_a (nA)$")
    plt.grid()
    plt.title(title)
    
    if save: plt.savefig("Imágenes/" + title + ".png", dpi=200)
    
def regLineal(x, y):
    
    aprox = np.polyfit(x, y, 1)
    a, b = aprox[0], aprox[1]
    
    sigma = np.sqrt(np.sum((y - a*x - b)**2)/(len(y)-2))
    dA = np.sqrt(len(y))*sigma/np.sqrt(len(y)*np.sum(x**2) - (np.sum(x)**2))
    dB = dA*np.sqrt(np.sum(x**2)/len(y))
    c = (np.sqrt(np.sum((x-np.mean(x))**2))*np.sqrt((np.sum((y-np.mean(y))**2))))
    r = np.sum((x-np.mean(x))*(y-np.mean(y)))/c
    
    return a, b, dA, dB, r**2

def UminLin(minimos, figs=(9, 6), save=False, title="Gráfica de mínimos"):
    
    n = np.array([i+1 for i in range(len(minimos))])
    nt = np.linspace(1, max(n), 1000)
    a, b, dA, dB, r2 = regLineal(n, minimos)
    aprox = np.array([a, b])
    
    plt.figure(figsize=figs)
    plt.scatter(n, minimos)
    plt.plot(nt, np.polyval(aprox, nt), "--r")
    plt.grid()
    plt.text(max(n)-1, max(minimos)- 2, "y = " + str(round(a, 2)) + "x + " +
             str(round(b, 2)) + "\n" + "$R^2 = $" + str(round(r2, 6))) 
    plt.xlabel("$n$")
    plt.ylabel("$U_{min} (V)$")
    plt.title(title)
    
    if save: plt.savefig("Imágenes/" + title + ".png", dpi=200)
    
    return a, b, dA, dB, r2

def deltaMin(minimos):
    
    DeltaMin = []
    
    for i in range(1, len(minimos)):
        
        DeltaMin.append(minimos[i] - minimos[i-1])
        
    return DeltaMin
    
def mediaMins(minimos):
    
    DeltaMin = deltaMin(minimos)
        
    return np.mean(DeltaMin), np.std(DeltaMin)

def recorridoLibre(minimos, save=False, title="Recorrido libre"):
    
    DeltaMin = deltaMin(minimos)
    x = np.array([2*i - 1 for i in range(1, len(DeltaMin)+1)])
    y = DeltaMin
    aprox = np.polyfit(x, y, 1)
    
    plt.figure()
    plt.grid()
    plt.title(title)
    plt.xlabel("$2n - 1$")
    plt.ylabel("$\\Delta U_{1, min} (eV)$")
    plt.scatter(x, y)
    plt.plot(x, np.polyval(aprox, x), "--")
    
    if save: plt.savefig("Imágenes/" + title + ".png", dpi=200)
    
    return regLineal(x, y)
    
path = "Datos Neon/Feo.txt"
U1, Ia, minimos = dataFH(path)
plotIa(path, save=True, title="Ne Feo")
print(mediaMins(minimos))
print(UminLin(minimos))
print(recorridoLibre(minimos))