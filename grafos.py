WHITE = "white"
BLACK = "black"
GRAY = "gray"
from copy import copy
import math


class Graph:
    def _buildAdjMatrix(self):
        self.adjMat = [[0 for v in range(len(self.vertexes))] for v in range(len(self.vertexes))]  # V2
        for relation in self.relations:
            row, col = self.encoder[relation[0]], self.encoder[relation[1]]
            self.adjMat[row][col] = 1

    def _buildEncoding(self):
        self.encoder, self.decoder = {}, {}
        index = 0
        for v in self.vertexes:
            self.encoder[v] = index
            self.decoder[index] = v
            index = index + 1

    def _buildAdjList(self):
        self.adjList = {}
        for v in self.vertexes:
            self.adjList[v] = []
        for relation in self.relations:
            self.adjList[relation[0]].append(relation[1])

    def _buildRelation(self, e):
        if self.directed:
            self.relations = e
        else:
            self.relations = set()
            for el in e:
                self.relations.add(el)
                self.relations.add((el[1], el[0]))

    def __init__(self, v, e, directed=True, view=True):
        self.directed = directed
        self.view = view  # Si True estoy usando listas de adyacencias / False, debo usar la matriz de adyacencias
        self.vertexes = v
        self._buildRelation(e)
        self._buildAdjList()
        self._buildEncoding()
        self._buildAdjMatrix()

    def getAdjMatrix(self):
        return self.adjMat

    def getAdjList(self):
        return self.adjList

    def _buildVProps(self, source=None):
        self.v_props = {}
        for v in self.vertexes:
            self.v_props[v] = {
                'color': WHITE,
                'distance': math.inf,
                'parent': None
            }
        if source is not None:
            self.v_props[source] = {
                'color': GRAY,
                'distance': 0,
                'parent': None
            }

    def _getNeighborsAdjList(self, vertex):
        return self.adjList[vertex]

    def _getNeighborsMatAdj(self, vertex):
        fila = self.adjMat[vertex].copy()
        vecinos = []
        while 1 in fila:
            vecinos.append(fila.index(1))
            fila[fila.index(1)] = 0
        return vecinos  # Completar esta función para el laboratorio.

    def getNeighbors(self, vertex):
        if self.view:
            return self._getNeighborsAdjList(vertex)
        # Retornar los vecinos de forma correcta para la representación de la matriz de adyacencias
        return self._getNeighborsMatAdj(vertex)

    def bfs(self, source):
        self._buildVProps(source)
        queue = [source]
        while len(queue) > 0:
            u = queue.pop(0)
            for neighbor in self.getNeighbors(u):
                if self.v_props[neighbor]['color'] == WHITE:
                    self.v_props[neighbor]['color'] = GRAY
                    self.v_props[neighbor]['distance'] = self.v_props[u]['distance'] + 1
                    self.v_props[neighbor]['parent'] = u
                    queue.append(neighbor)
            self.v_props[u]['color'] = BLACK
        return self.v_props

    def dfs(self):
        self._buildVProps()
        time = 0
        for v in self.vertexes:
            if self.v_props[v]['color'] == WHITE:
                time = self.dfs_visit(v, time)
        return self.v_props

    def dfs_visit(self, vertex, time):
        time = time + 1
        self.v_props[vertex]['distance'] = time
        self.v_props[vertex]['color'] = GRAY
        for neighbor in self.getNeighbors(vertex):
            if self.v_props[neighbor]['color'] == WHITE:
                self.v_props[neighbor]['parent'] = vertex
                time = self.dfs_visit(neighbor, time)
        self.v_props[vertex]['color'] = BLACK
        time = time + 1
        self.v_props[vertex]['final'] = time
        return time


def printZonas(v_props, grafo):
    for v in v_props.keys():
        if grafo.v_props[v]['color'] == BLACK:
            v_props[v]['path'] = '-->'.join(map(str, getPath(v, v_props)))
            print(str(v), '-->', v_props[v])


def printVProps(v_props, grafo=None):
    for v in v_props.keys():
        v_props[v]['path'] = '-->'.join(map(str, getPath(v, v_props)))
        print(str(v), '-->', v_props[v])


def printAdjMatrix(graph):
    print("===================== ADJ Matrix =======================")
    adjMat = graph.getAdjMatrix()
    for row in adjMat:
        print(' '.join(list(map(str, row))))


def getPath(vertex, v_props):
    path = [vertex]
    current = vertex
    while (v_props[current]['parent'] is not None):
        path.insert(0, v_props[current]['parent'])
        current = v_props[current]['parent']
    return path


def printAdjList(graph):
    print("===================== ADJ List =======================")
    adjList = graph.getAdjList()
    for v in adjList.keys():
        print(str(v), list(map(str, adjList[v])))


def checkIfConnectedComponent(result_bfs):
    amount_black = 0
    for e in result_bfs.keys():
        if result_bfs[e]['color'] == BLACK:
            amount_black += 1
    return amount_black > 1


def main():
    print("********************* Caso de negociación  *********************")
    # Como equipo hemos escogido usar un servicio de mensajería de la empresa servientrega.
    # Para este ejemplo se van usar dos distrbuidoras, una ubicada en Medellin y la ota ubicada en Bogotá.
    # Para este caso se van a usar BFS y DFS  para determinar dependiendo de las zonas que hayan realizado pedido en el mes, que distribuidora va a ser el envio y cual es la ruta más cercana para cubrir todos sus destinos
    # Sea:

    CGP = "Centro de gestión en el Poblado"
    Pl = "Popular"
    Scz = "Santa Cruz"
    M = "Manrique"
    Ar = "Aranjuez"
    Ct = "Castilla"
    Oc = "12 de Octubre"
    Rb = "Roblado"
    Vh = "Villa hermosa"
    Ba = "Buenos aires"
    Lcm = "La candelaria"
    L = "Laureles"
    La = "La américa"
    Sj = "San javier"
    G = "Guayabal"
    B = "Belén"

    CGK = "Centro de gestión de Servientrega en Kennedy"
    Mr = "Los Martires"
    F = "Fontibon"
    Sf = "Santa fé"
    An = "Antonio Nariño"
    Pa = "Puente Aranda"
    T = "Teusaquillo"
    Tj = "Tunjuelito"
    Cb = "Ciudad Bolivar"
    Us = "Usme"
    Sc = "San Cristobal"
    LCB = "La candelaria"
    Ch = "Chapinero"
    U = "Usaquen"
    Sb = "Suba"
    Bu = "Barrios Unidos"
    Ru = "Rafael Uribe"
    Bs = "Bosa"
    E = "Engativa"

    zonas = {CGK, Mr, F, Sf, An, Pa, T, Cb, Us, Sc, LCB, Ch, Tj, U, Sb, Bu, Ru, Bs, CGP, Pl, Scz, M, Ar, Ct, Oc, Rb, Vh,
             Ba, Lcm, L, La, Sj, G, B, E}
    relaciones = {(Pl, Scz), (Pl, M), (Scz, Pl), (Scz, Ct), (Scz, Ar), (M, Pl), (M, Ar), (M, Vh), (Ar, Scz), (Ar, Lcm),
                  (Ar, Ct), (Ar, M), (Ct, Scz), (Ct, Ar), (Ct, Oc), (Ct, Rb), (Oc, Ct), (Oc, Rb), (Rb, Oc), (Rb, L),
                  (Rb, Sj), (Rb, Ct), (Vh, M), (Vh, Lcm), (Vh, Ba), (Ba, Vh), (Ba, Lcm), (Ba, CGP), (Lcm, Ar),
                  (Lcm, CGP), (Lcm, Ba), (Lcm, Ba), (L, Lcm), (L, Rb), (L, B), (L, La), (La, Sj), (La, L), (Sj, La),
                  (Sj, Rb), (CGP, G), (CGP, Ba), (CGP, Lcm), (G, CGP), (G, B), (B, G), (B, L), (CGK, Bs), (CGK, Cb),
                  (CGK, Tj), (CGK, Pa), (CGK, F), (Mr, An), (Mr, Sf), (Mr, T), (Mr, Pa), (F, E), (F, CGK), (Sf, LCB),
                  (Sf, Ch), (Sf, Sc), (Sf, Mr), (An, Mr), (An, Ru), (An, Sc), (An, Tj), (An, Pa), (Pa, Mr), (Pa, Tj),
                  (Pa, An), (Pa, CGK), (Pa, T), (T, Ch), (T, Bu), (T, E), (T, Mr), (T, F), (T, Pa), (Tj, An), (Tj, Cb),
                  (Tj, Bs), (Tj, Ru), (Tj, CGK), (Tj, Pa), (Cb, Tj), (Cb, Bs), (Cb, Us), (Us, Sc), (Us, Ru), (Us, Cb),
                  (Sc, Sf), (Sc, Us), (Sc, An), (Sc, Ru), (LCB, Sf), (Ch, U), (Ch, Bu), (Ch, T), (Ch, Sf), (U, Ch),
                  (U, Sb), (Sb, U), (Sb, E), (Sb, Bu), (Bu, Ch), (Bu, T), (Bu, E), (Bu, Sb), (Ru, U), (Ru, Tj),
                  (Ru, An), (Ru, Sc), (Bs, CGK), (Bs, Cb), (Bs, Tj), (E, F), (E, Sb), (E, T), (E, Bu)}

    rutas = Graph(zonas, relaciones, False)
    printAdjList(rutas)
    printAdjMatrix(rutas)
    # Se utiliza el BFS con S = CGP y luego S = CGK para determinar que zonas cubre en cada ciudad
    print("====================== Zonas a las cuales distribuye servientrega en Medellin ===================")
    printZonas(rutas.bfs(CGP), rutas)
    print("====================== Zonas a las cuales distribuye servientrega en Bogotá ===================")
    printZonas(rutas.bfs(CGK), rutas)
    print("====================== Busca por profundidad primero ===================")
    result = rutas.dfs()
    printVProps(result)


main()
