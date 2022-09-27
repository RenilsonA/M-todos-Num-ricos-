import math as m
import matplotlib.pyplot as plt

alternativa = ["  A) ", "  B) ", "  C) ", "  D) ", "  E) "]
pi = 3.141592653589793
e = 2.718281828459045

def converteDecimal(num):
    espacos = 0
    bitsExpoente = 11
    bitsMantissa = 0
    sinal = 1
    expoente = 0
    fracao = 0
    for i in range(len(num)):
        if num[i] == ' ':
            espacos = espacos + 1
        elif espacos == 2:
            bitsMantissa = bitsMantissa + 1
            fracao = fracao + (ord(num[i]) - 48)*pow(1/2, bitsMantissa)
        elif espacos == 1:
            bitsExpoente = bitsExpoente - 1
            expoente = expoente + (ord(num[i]) - 48)*pow(2, bitsExpoente)
        elif espacos == 0:
            if num[i] == '1':
                sinal = -1
    return sinal * pow(2, expoente - 1023) * (fracao + 1)

def fatorial(num):
    x = 1
    while num > 1:
        x = num * x
        num = num - 1
    return x

def positivo(num):
    if num < 0:
        num = -num
    return num
    
def ordemdoNum(num):
    div = 1
    num = positivo(num)
    while num/div >= 1:
        div = div * 10
    return div

def truncar(num, casas):
    x = pow(10, casas)
    div = ordemdoNum(num)
    x = x/div
    return int(num*x)/x

def arredondamento(num, casas):
    x = pow(10, casas)
    num = truncar(num, casas + 1)
    adicionar = 0.0
    div = ordemdoNum(num)
    y = 1
    x = x/div
    if num < 0:
        y = -1
    if ((num*x)%y)*10 >= 5.0 or ((num*x)%y)*10 <= -5.0:
        adicionar = y/x
    return int(num*x)/x + adicionar

def Raizdelta(a, b, c):
    return pow(b * b - 4 * a * c, 0.5)

def raiz2grau(a, b, c):
    D = Raizdelta(a, b, c)
    x1 = (-b + D) / (2 * a)
    x2 = (-b - D) / (2 * a)
    return x1, x2

def erroAbsoluto(x1, x2):
    return (positivo(x1 - x2))

def erroRelativo(EA, x):
    return EA/positivo(x)

def somatorio4C(i, valor):
    x1 = 0
    while i > 0:
        x1 = x1 + (pow(-1, i)*pow(valor, i)/fatorial(i))
        i = i - 1
    return x1

#entrada: 5x2^(2x)sin(2x)ln(3x)
#i[0] = mult da funcao          i[1] = exp de X             i[2] = base da exp por X       i[3] = mult do exp vezes X     
#i[4] = func trigonometrica     i[5] = mult da func trig    i[6] = há ln                   i[7] = mult interno do ln    
#i[0] = 5, i[1] = 1, i[2] = 2, i[3] = 2, i[4] = True = seno, i[5] = 2, i[6] = True, i[7] = 3 
def funcao(equacao, x):
    f = 0
    for i in equacao:
        func = i[0]*pow(x, i[1])
        if func == 0:
            func = 1
        if len(i) >= 4:
            func = func*pow(i[2], i[3] * x)
        if len(i) >= 6:
            if i[4] == 1:
                func = func*m.sin(i[5]*x)
            elif i[4] == 2:
                func = func*m.cos(i[5]*x)
        if len(i) >= 8:
            if i[6]:
                func = func*m.log(i[7]*x, e)
        f = f + func
    return f

#funciona apenas para algumas equações simples e especialmente para as alternativas da questão 9
def derivada(equacao, x): 
    f = 0
    for i in equacao:
        func = 1
        tam = len(i)
        if x != 0:
            func = i[0]*i[1]*pow(x, i[1] - 1)
        if func == 0 and tam > 2:
            func = 1
        if tam >= 4 and i[3] != 0:
            func = - func*i[3]*m.log(i[2], e)*pow(i[2], i[3]*x)
        if tam >= 6:
            if i[1] != 0:
                if i[4] == 1:
                    func = func*(m.sin(i[5]*x) + i[5]*x*m.cos(i[5]*x))
                elif i[4] == 2:
                    func = func*(m.cos(i[5]*x) - i[5]*x*m.sin(i[5]*x))
            elif i[4] == 1:
                func = func*m.cos(i[5]*x)*i[5]*i[0]
            elif i[4] == 2:
                func = - func*m.sin(i[5]*x)*i[5]*i[0]
        if tam >= 8:
            if i[6]:
                func = func*(1/(i[7]*x))*i[7]
        f = f + func
    return f

def biseccao(x1, x2, equacao, tolerancia, iteracoes, pn):
    condicao = True
    x3 = 0
    valorPN, valorYPN = [], []
    if funcao(equacao, x1) * funcao(equacao, x2) > 0:
        print("valores da função inválidos")
        return
    while condicao and iteracoes > 0:
        x3 = (x1 + x2)/2
        y1 = funcao(equacao, x1)
        y3 = funcao(equacao, x3)
        if y1 * y3 < 0:
            x2 = x3
        else:
            x1 = x3
        valorPN.append(x3)
        valorYPN.append(funcao(equacao, x3))
        condicao = abs(funcao(equacao, x3)) > tolerancia
        iteracoes = iteracoes - 1
    if pn:
        return valorPN, valorYPN
    return x3

def newton(x, equacao, tolerancia, iteracoes, pn):
    condicao = True
    listaPN, listaYPN = [], []
    while condicao and iteracoes > 0:
        if derivada(equacao, x) == 0:
            print("Erro, dividindo por zero.")
            return
        x1 = x - funcao(equacao, x) / derivada(equacao, x)
        listaPN.append(x1)
        listaYPN.append(funcao(equacao, x1))
        if positivo(x1 - x) < tolerancia:
            if pn:
                return listaPN, listaYPN
            return x1
        x = x1
        iteracoes = iteracoes - 1
    print("Metodo de Newton falhou, iterações insuficientes.")
    return

def graficoPN(x1, x2, equacao, tolerancia, iteracoes, tipo):
    pn, ypn = [], []
    if tipo == 1:
        pn, ypn = biseccao(x1, x2, equacao, tolerancia, iteracoes, True)
    else:
        pn, ypn = newton((x1 + x2)/2, equacao, tolerancia, iteracoes, True)
    if pn:
        plt.scatter(pn, ypn, color='gray', label = 'Aproximações')
        plt.scatter(pn[-1], ypn[-1], color='red', label = 'Raiz')

def grafico(x1, x2, equacao, delta):
    valoresx = []
    valoresy = []
    if x1 > x2:
        x3 = x1
        x1 = x2
        x2 = x3
    while x1 < x2:
        valoresx.append(x1)
        valoresy.append(funcao(equacao, x1))
        x1 = x1 + delta
    plt.plot(valoresx, valoresy, color='green', label = 'F(x)')
    plt.axvline(0, color='black')
    plt.axhline(0, color='black')
    plt.legend()                  
    plt.show()

def pontoFixo(A, p0, tolerancia, N):
    p = 0
    i = N
    while i > 0:
        p = 2*p0 - A*p0*p0
        if positivo(p - p0) < tolerancia:
            return p
        p0 = p
        i = i - 1
    print("Falhou")
    return -1

def questao1():
    valores = ["0 10000001010 1001001100000000000000001000000000001000000000000000",
               "1 10000001010 1001001100000000001000000000000000000000000000000000",
               "0 01111011111 0101001100000000000000000000000000000000000001000000",
               "0 01011111011 0101001100000000000000000000000000000000001000000001"]
    altern = 0
    for i in valores:
        valor = converteDecimal(i)
        print(alternativa[altern], valor)
        altern = altern + 1

def questao2():
    print("Questão 2:")
    dadosP  = [pi, pow(2, 0.5), pow(e, 10), fatorial(9)]
    dadosP2 = [3, 1.414, 22000, (pow(pi*18, 0.5)*pow((9/e), 9))]
    for i in range(4):
        EA = erroAbsoluto(dadosP[i], dadosP2[i])
        ER = erroRelativo(EA, dadosP[i])
        print(alternativa[i] + "Erro absoluto de", EA, "e erro relativo de", ER)

def questao3():
    print("Questão 3:")
    dadosP  = [150, 900, 1500, 90]
    for i in range(4):  
        # (|p - p*|/|p|) <= 10^-3 >> |p - p*| <= p*10^-3 >> - p - p*10^-3 <= - p* <= - p + p*10^-3 >>
        # p + p*10^-3 <= p* <= p - p*10^-3 
        LInfP2 = dadosP[i] + dadosP[i]*pow(10, -3)
        LSnfP2 = dadosP[i] - dadosP[i]*pow(10, -3)
        print(alternativa[i] + "O intervalo de p* é de:", LInfP2, "<= p* <=", LSnfP2)

def questao4():
    print("Questão 4:")
    x = -10*pi + 6*e - 0.327
    xtrunc3 = -10*truncar(pi, 3) + 6*truncar(e, 3) - 0.327
    xtrunc4 = -10*truncar(pi, 4) + 6*truncar(e, 4) - 0.327
    xtrunc4 = truncar(xtrunc4, 4)
    xarred3 = -10*arredondamento(pi, 3) + 6*arredondamento(e, 3) - 0.327
    xarred4 = -10*arredondamento(pi, 4) + 6*arredondamento(e, 4) - 0.327
    xarred4 = arredondamento(xarred4, 4)
    EA = erroAbsoluto(x, xtrunc3)
    ER = erroRelativo(EA, x)
    print(alternativa[0] + "Truncamento de 3 casas: erro absoluto de", EA, "e erro relativo de", ER)
    EA = erroAbsoluto(x, xtrunc4)
    ER = erroRelativo(EA, x)
    print("   Truncamento de 4 casas: erro absoluto de", EA, "e erro relativo de", ER)
    EA = erroAbsoluto(x, xarred3)
    ER = erroRelativo(EA, x)
    print("   Arredondamento de 3 casas: erro absoluto de", EA, "e erro relativo de", ER)
    EA = erroAbsoluto(x, xarred4)
    ER = erroRelativo(EA, x)
    print("   Arredondamento de 3 casas: erro absoluto de", EA, "e erro relativo de", ER)

    a, b, c = 1/3, -123/4, 1/6
    x1, x2 = raiz2grau(a, b, c)
    At3, At4 = truncar(a, 3), truncar(a, 4)
    Bt3, Bt4 = truncar(b, 3), truncar(b, 4)
    Ct3, Ct4 = truncar(c, 3), truncar(c, 4)
    Aa3, Aa4 = arredondamento(a, 3), arredondamento(a, 4)
    Ba3, Ba4 = arredondamento(b, 3), arredondamento(b, 4)
    Ca3, Ca4 = arredondamento(c, 3), arredondamento(c, 4)
    x1t3, x2t3 = raiz2grau(At3, Bt3, Ct3)
    x1t4, x2t4 = raiz2grau(At4, Bt4, Ct4)
    x1a3, x2a3 = raiz2grau(Aa3, Ba3, Ca3)
    x1a4, x2a4 = raiz2grau(Aa4, Ba4, Ca4)
    x1trunc3, x2trunc3 = truncar(x1t3, 3), truncar(x2t3, 3)
    x1trunc4, x2trunc4 = truncar(x1t4, 4), truncar(x2t4, 4)
    x1arred3, x2arred3 = arredondamento(x1a3, 3), arredondamento(x2a3, 3)
    x1arred4, x2arred4 = arredondamento(x1a4, 4), arredondamento(x2a4, 4)
    EA = erroAbsoluto(x1, x1trunc3)
    ER = erroRelativo(EA, x1)
    EA2 = erroAbsoluto(x2, x2trunc3)
    ER2 = erroRelativo(EA2, x2)
    print(alternativa[1] + "Truncamento de 3 casas: erro absoluto de X1", EA, "e de X2", EA2, "e erro relativo de X1", ER, "e de X2", ER2)
    EA = erroAbsoluto(x1, x1trunc4)
    ER = erroRelativo(EA, x1)
    EA2 = erroAbsoluto(x2, x2trunc4)
    ER2 = erroRelativo(EA2, x2)
    print("     Truncamento de 4 casas: erro absoluto de X1", EA, "e de X2", EA2, "e erro relativo de X1", ER, "e de X2", ER2)
    EA = erroAbsoluto(x1, x1arred3)
    ER = erroRelativo(EA, x1)
    EA2 = erroAbsoluto(x2, x2arred3)
    ER2 = erroRelativo(EA2, x2)
    print("     Arredondamento de 3 casas: erro absoluto de X1", EA, "e de X2", EA2, "e erro relativo de X1", ER, "e de X2", ER2)
    EA = erroAbsoluto(x1, x1arred4)
    ER = erroRelativo(EA, x1)
    EA2 = erroAbsoluto(x2, x2arred4)
    ER2 = erroRelativo(EA2, x2)
    print("     Arredondamento de 3 casas: erro absoluto de X1", EA, "e de X2", EA2, "e erro relativo de X1", ER, "e de X2", ER2)

    x = pow(e, -5)
    x1 = somatorio4C(9, 5)
    xtrunc3 = truncar(x1, 3)
    xtrunc4 = truncar(x1, 4)
    xarred3 = arredondamento(x1, 3)
    xarred4 = arredondamento(x1, 4)
    print(xtrunc3, xtrunc4, xarred3, xarred4)
    EA = erroAbsoluto(x, xtrunc3)
    ER = erroRelativo(EA, x)
    print(alternativa[2] + "Truncamento de 3 casas: erro absoluto de X1", EA, "e erro relativo de X1", ER)
    EA = erroAbsoluto(x, xtrunc4)
    ER = erroRelativo(EA, x)
    print("     Truncamento de 4 casas: erro absoluto de X1", EA, "e erro relativo de X1", ER)
    EA = erroAbsoluto(x, xarred3)
    ER = erroRelativo(EA, x)
    print("     Arredondamento de 3 casas: erro absoluto de X1", EA, "e erro relativo de X1", ER)
    EA = erroAbsoluto(x, xarred4)
    ER = erroRelativo(EA, x)
    print("     Arredondamento de 3 casas: erro absoluto de X1", EA, "e erro relativo de X1", ER)

def questao5():
    print("As funções estão no cabeçalho, com nomes de função (deve-se digitar nos vetores acima os valores")
    print("respectivos dos termos multiplicadores e seus graus, em constantes, colocar como 0), derivada,")
    print("biseccao(último parametro pode por false por enquanto, pois ele define o retorno de pn's da questão")
    print("7) e newton (último parâmetro como False também).")

def questao6():
    print("A função está no cabeçalho e é chamada de grafico, onde ela calcula cada ponto de y dado um intervalo de x")

def questao7():
    print("As funções de bisseccao e newton já fazem isso, basta colocar o termo pn como True, que eles retornam os valores de \
           pn e Y(pn).")

def questao8():
    print("As funções de graficoPN do cabeçalho resolve, onde deve-se passar o intervado [x1, x2], equação, a tolerância e \
           o número máximo de iterações. Essa função já pega os valores de pn e Y(pn), logo não precisa definir nada.")

def questao9():
    #Comentários da linha 70 para entender melhor.
    equacao = [((1, 1), (-1, 0, 2, -1)), 
               ((1, 1), (1, 0), (-2, 0, 1, 0, 1, pi)), 
               ((2, 1, 1, 0, 2, 2), (-1, 2), (-2, 1), (-1, 0)), 
               ((1, 0, 1, 0, 0, 0, True, 1), (-1, 0, 2, 1), (1, 2))]
    intervalo = [(0, 1), (0, 0.5), (-3, -2), (3, 5)]
    tolerancia = pow(10, -4)
    iteracoes = 20
    func = ["Bisecção", "Newton"]
    #Tipo 1 = biseccao (muitos pontos, mas estável), Tipo 2 = Newton (menos pontos)
    tipo = 2
    delta = 0.001
    for i in range(len(equacao)):
        valorX3 = biseccao(intervalo[i][0], intervalo[i][1], equacao[i], tolerancia, iteracoes, False)
        valorX2 = newton((intervalo[i][0] + intervalo[i][1]) / 2, equacao[i], tolerancia, iteracoes, False)
        valorY3 = funcao(equacao[i], valorX3)
        valorY2 = funcao(equacao[i], valorX2)
        pn, ypn = newton((intervalo[i][0] + intervalo[i][1]) / 2, equacao[i], tolerancia, iteracoes, True)
        
        print(alternativa[i])
        print("     Valor de pn e F(pn) em bisecção:", valorX3, "e", valorY3)
        print("     Valor de pn e F(pn) na função Newton:", valorX2, "e", valorY2)
        print("     Valores de pn e Y(pn), com base na função ", func[tipo - 1], pn, ypn)

        graficoPN(intervalo[i][0], intervalo[i][1], equacao[i], tolerancia, iteracoes, tipo)
        grafico(intervalo[i][0], intervalo[i][1], equacao[i], delta)

def questao10():
    print(alternativa[0], "Se g(x) = 2x-Ax², g(p) = p = 2p-Ap², e assumindo que A > 0")
    print("     Encontrando as raizes da equação algebricamente, temos que:")
    print("     2p - Ap² - p = 0")
    print("     p(1 - Ap) = 0")
    print("     Extraindo as raízes")
    print("     p' = 0")
    print("     1 - Ap = 0")
    print("     p = 1/A")
    print("     Como A > 0, logo p converge.")
    print("     Logo a iteração é válida.")

    A = 2
    p = 0.99
    p = pontoFixo(A, p, 0.0001, 100)
    print(alternativa[1], p)

    print(alternativa[2], "Igualando a equação dada a zero, podemos encontrar as raizes da equação algebricamente:")
    print("     2p - Ap² = 0")
    print("     p(2 - Ap) = 0")
    print("     Extraindo as raízes")
    print("     p' = 0")
    print("     2 - Ap = 0")
    print("     p'' = 2/A")
    print("     Logo, temos que no invervalo [0, 2/A] existe uma raiz válida.")
    print("     Então, podemos afirmar que o intervalo varia de 0 <= g(x) <= 2/A")

    
questao10()
