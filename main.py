from tkinter import *
from tkinter.ttk import Combobox, Style
from tabulate import tabulate
import numpy as np
import matplotlib.pyplot as plt

"""
Пример параметров (неявная разностная схема):
    I = 50
    K = 100
    T = 200 c
    R = 8 см
    betta = 0.02 1/см
"""

def clicked():
    #сбор параметров
    I = int(txt_I.get())
    K = int(txt_K.get())
    T = float(txt_T.get())
    R = float(txt_R.get())
    B = float(txt_B.get())
    flag = selected.get()
    shem(I, K, R, T, B, flag)
    #выбор схемы 

          
def Ir1(r, a_0):
    return 250 * np.exp(-(r / a_0) ** 2)


def shem(I, K, R, T, B, flag):
    w = calculate(I, K, T, R, B)
    h_t = T / (K - 1)
    h_r = R / (I - 1)
    drawNormal(h_t, h_r, K, w, T, R, I, flag)
    errorNormal(I, K, T, R, B)
    
   
def tridiagonal_matrix_algorithm(a, b):
    n = len(a)
    x = [0 for k in range(0, n)]
    # Прямой ход
    v = [0 for k in range(0, n)]
    u = [0 for k in range(0, n)]
    # для 0-й строки
    v[0] = -a[0][1] / a[0][0]
    u[0] = b[0] / a[0][0]
    for i in range(1, n - 1):  # заполняем за исключением 0-й и (n-1)-й строк матрицы
        v[i] = -a[i][i + 1] / (a[i][i] + a[i][i - 1] * v[i - 1])
        u[i] = (-a[i][i - 1] * u[i - 1] + b[i]) / (a[i][i] + a[i][i - 1] * v[i - 1])
    # для (n-1)-й строки
    v[n - 1] = 0
    u[n - 1] = (-a[n - 1][n - 2] * u[n - 2] + b[n - 1]) / (a[n - 1][n - 1] + a[n - 1][n - 2] * v[n - 2])
    # Обратный ход
    x[n - 1] = u[n - 1]
    for i in range(n - 1, 0, -1):
        x[i - 1] = v[i - 1] * x[i] + u[i - 1]

    return x


def calculate(I, K, T, R, B):
    k = 0.065
    c = 1.35
    alpha = 0.001
    L = 1
    w_0 = 0
    a_0 = R / 8
    w = np.zeros((I, K))
    b1 = np.zeros(I)
    b = np.zeros(I)

    h_t = T / (K - 1)
    h_r = R / (I - 1)
    gamma = k * h_t / h_r ** 2
    coef = 2 * alpha * h_t / L
    P = np.zeros((I, I), np.float64)

    P[0][0] = c + 4 * gamma + coef
    P[0][1] = -4 * gamma
    P[I - 1][I - 1] = c + gamma + coef - h_r / (2 * R)
    P[I - 1][I - 2] = h_r / (2 * R) - gamma
    P[I - 2][I - 2] = c + gamma + coef - h_r / (2 * (R - h_r))
    P[I - 2][I - 3] = h_r / (2 * (R - h_r)) - gamma


    for i in range(1, I - 2):
        P[i][i - 1] = gamma * (h_r / (2 * i * h_r) - 1)
        P[i][i] = c + 2 * gamma + coef
        P[i][i + 1] = -gamma * (h_r / (2 * i * h_r) + 1)

    for i in range(0, I):
        w[i][0] = w_0
        b1[i] = h_t * B * Ir1(i * h_r, a_0)

    for k in range(1, K):
        b = [b[k - 1] * c for b in w] + b1
        ans = tridiagonal_matrix_algorithm(P, b)
        for i in range(0, I - 1):
            w[i][k] = ans[i]
        w[I - 1][k] = w[I - 2][k]

    return w


#погрешности
def errorNormal(I, K, T, R, B):
    epsilon = calculate(I, K, T, R, B)[int(I/2)][int(K/2)] - calculate(I, K*2-1, T, R, B)[int(I/2)][K]
    epsilon2 = calculate(I, K*2-1, T, R, B)[int(I/2)][K] - calculate(I, K*4-3, T, R, B)[int(I/2)][K*2]
    value_list = list()
    value_list.append([I, K, round(abs(epsilon), 5),round(abs(epsilon2), 5), round(abs(epsilon/epsilon2), 5)])
    column_list = ["I", "K", "∆_w_ht,hr", "∆_w_ht/2,hr", "delta_ht,hr"]
    print(tabulate(value_list, column_list, tablefmt="grid"))
      
        
def drawNormal(h_t, h_r, K, w, T, R, I, flag):
    if flag == 1:
        fig1 = plt.figure()
        ax = fig1.add_subplot()
        ax.plot([i * h_t for i in range(0, K)], w[int(1 / h_r) - 1][:], label='r = 1 см')
        ax.plot([i * h_t for i in range(0, K)], w[int(2 / h_r) - 1][:], label='r = 2 см')
        ax.plot([i * h_t for i in range(0, K)], w[int(3 / h_r) - 1][:], label='r = 3 см')
        ax.plot([i * h_t for i in range(0, K)], w[int(5 / h_r) - 1][:], label='r = 5 см')
        ax.plot([i * h_t for i in range(0, K)], w[int(8 / h_r) - 1][:], label='r = 8 см')
        ax.grid()
        ax.set_xlabel('T, c')
        ax.set_ylabel('w, град')
        ax.legend()
    elif flag == 0:
        fig2 = plt.figure()
        ax = fig2.add_subplot()
        ax.plot([i * h_r for i in range(0, I)], [x[int(15 / h_t) - 1] for x in w], label='t = 15с')
        ax.plot([i * h_r for i in range(0, I)], [x[int(50 / h_t) - 1] for x in w], label='t = 50с')
        ax.plot([i * h_r for i in range(0, I)], [x[int(100 / h_t) - 1] for x in w], label='t = 100с')
        ax.plot([i * h_r for i in range(0, I)], [x[int(150 / h_t) - 1] for x in w], label='t = 150с')
        ax.plot([i * h_r for i in range(0, I)], [x[int(200 / h_t) - 1] for x in w], label='t = 200с')
        ax.grid()
        ax.set_xlabel('r, cм')
        ax.set_ylabel('w, град')
        ax.legend()
    plt.show()


window = Tk()
window.title("Курсовая ЧММФ")
"""
window.style = Style(window)
window.style.theme_use("xpnative")
"""
frame = Frame(window, borderwidth=1)
frame.pack(fill=X)
lbl = Label(frame, text="Численное моделирование динамического теплового поля оптического элемента, освещаяемого лазером.", font="Arial 14", padx=10, pady=10)
lbl.pack(side = TOP)
lbl2 = Label(frame, text="Неявная разностная схема.", font="Arial 13", padx=10)
lbl2.pack(side = LEFT)

frameGrid = Frame(window)
frameGrid.pack(fill=X)
lbl3 = Label(frameGrid, text="Параметры сетки:", font="Arial 13")
lbl3.pack(side = LEFT, padx=10, pady=10)

frame1 = Frame(window)
frame1.pack(fill=X)
lbl4 = Label(frame1, text="Число строк I: ", font="Arial 13", padx=10, pady=5)
lbl4.pack(side = LEFT)
txt_I = Entry(frame1,width=20)
txt_I.place(x = 160, y = 7)

frame2 = Frame(window)
frame2.pack(fill=X)
lbl5 = Label(frame2, text="Число столбцов K: ", font="Arial 13", padx=10, pady = 5)
lbl5.pack(side = LEFT)
txt_K = Entry(frame2,width=20)
txt_K.place(x = 160, y = 7)

frame3 = Frame(window)
frame3.pack(fill=X)
lbl6 = Label(frame3, text="Параметры системы:", font="Arial 13", padx=10, pady = 10)
lbl6.pack(side = LEFT)

frame4 = Frame(window)
frame4.pack(fill=X)
lbl7 = Label(frame4, text="Время T [c]:",font="Arial 13", padx=10, pady = 5)
lbl7.pack(side = LEFT)
txt_T = Entry(frame4,width=20)
txt_T.place(x = 310, y = 7)

frame5 = Frame(window)
frame5.pack(fill=X)
lbl8 = Label(frame5, text="Радиус оптического элемента R [см]:",font="Arial 13", padx=10, pady = 5)
lbl8.pack(side = LEFT)
txt_R = Entry(frame5,width=20)
txt_R.place(x = 310, y = 7)

frame6 = Frame(window)
frame6.pack(fill=X)
lbl9 = Label(frame6, text="Коэффициент поглощения β [1/см]:",font="Arial 13", padx=10, pady = 5)
lbl9.pack(side = LEFT)
txt_B = Entry(frame6, width=20)
txt_B.place(x = 310, y = 7) 


frame7 = Frame(window)
frame7.pack(fill=X)
lbl20= Label(frame7, text="Построить график при фиксированном:",font="Arial 13", padx=10, pady = 5)
lbl20.pack(side = LEFT)
selected = IntVar()
selected.set(0)
rad1 = Radiobutton(window, text='I', value=0, variable=selected)
rad2 = Radiobutton(window, text='K', value=1, variable=selected) 
rad1.pack()
rad2.pack()


btn = Button(frame7, text="Посчитать", width=15, font="Arial 13", command=clicked)
btn.pack(side = RIGHT, padx=10, pady = 5, fill=Y)
window.mainloop()
