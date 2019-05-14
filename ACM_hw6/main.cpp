#include <limits>
#include <cstdint>
#include <iostream>
#include <vector>
#include <math.h>
#include <time.h>


/*
При реализации алгоритма Дейкстры редко используется фибоначчиева куча, дающая асимптотику O(E + V \log V),
из-за большой контанты. Реально же используются кучи, у которых 2, 4 или 8 потомков. Обычно 8-кучи оказываются
самыми эффективными из-за того, что запрашиваемые из памяти данные достаточно локальны.

Реализуйте структуру-хранилище очередь с приоритетами, интерфейс которой описан ниже.
*/

struct cmp : public std::unary_function<int, bool>
{
    explicit cmp(const int &a) : a(a) {}
    bool operator() (const int &arg)
    { return a == arg; }
    int a;
};

class DijkstraStorage {
    std::vector<int64_t> distances {};
    std::vector<int> heap {};
    size_t d = 0;
public:
    // Инициализирует структуру для графа на n_vertex вершинах
    explicit DijkstraStorage(size_t n_vertex, size_t d)
    : d(d) {
        distances.resize(n_vertex);
        distances.assign(n_vertex, std::numeric_limits<int64_t >::max());
    }

    // Возвращает текущее найденное расстояние до вершины. Если оно еще не найдено - возвращается
    // std::numeric_limits<int64_t>::max();
    int64_t GetDistanceToVertex(int vertex_index) const {
        return distances[vertex_index];
    }

    // Сообщает о том, что некоторый путь до вершины уже найден, но он еще не оптимален
    // Это означает, что vertex_index на данный момент распологается в очереди (на краю волнового фронта)
    bool IsActive(int vertex_index) const {
        return std::find(heap.begin(), heap.end(), vertex_index) != heap.end();
    }

    // Делает вершину активной и возвращает true, если вершина уже была активной возвращается false
    bool MakeVertexActive(int vertex_index) {
        if (IsActive(vertex_index))
            return false;

        heap.push_back(vertex_index);
        int i = heap.size() - 1;
        while (distances[heap[i]] < distances[heap[(int) floor(i/d)]]) {
            std::swap(heap[i], heap[(int) floor(i/d)]);
            i = (int) floor(i/d);
        }

        return true;
    }

    // Делает вершину не активной и возвращает true если вершина была активной, иначе возвращается false
    bool MakeVertexInactive(int vertex_index) {
        if (!IsActive(vertex_index))
            return false;

        distances[vertex_index] = -1;
        return true;
    }

    // Изменяет текущее расстояние до вершины, вершина может быть как активной так и не активной
    void SetDistance(int vertex_index, int64_t distance) {
        distances[vertex_index] = distance;
    }

    // Возвращает номер активной вершины, расстояние до которой минимально
    int GetMinIndex() const {
        return heap[0];
    }

    // Возвращает номер активной вершины, расстояние до которой минимально и делает ее неактивной
    int RemoveMin() {
        int min_index = heap[0];
        MakeVertexInactive(min_index);

        heap[0] = heap.back();
        heap.pop_back();

        int i = 0;
        while (d*i + d < heap.size()) {
            int64_t min_distance = std::numeric_limits<int64_t>::max();
            int best_index = 0;
            for (size_t j = d*i + 1; j <= d*i + d; j++) {
                if (distances[heap[j]] < min_distance) {
                    min_distance = distances[heap[j]];
                    best_index = heap[j];
                }
            }
            if (heap[i] < min_distance)
                break;
            std::swap(heap[i], heap[best_index]);
            i = best_index;
        }

        return min_index;
    }

    // true если не осталось активных вершин
    bool Empty() const {
        return heap.size() == 0;
    }

    void StorageDump() {
        std::cout <<"\nDistances: ";
        for(auto i = distances.begin(); i != distances.end(); i++)
            std::cout << *i << ' ';
        std::cout << std::endl;
        std::cout <<"Heap: ";
        for(auto i = heap.begin(); i != heap.end(); i++)
            std::cout << *i << ' ';
        std::cout << std::endl;
    }
};

class Graph {
public:
    int64_t **graph_matrix;
    size_t vertex_quan;
    Graph(size_t n_vertex);
    ~Graph();
    void AddEdge(size_t from, size_t to);
};

Graph::Graph(size_t n_vertex) {
    vertex_quan = n_vertex;
    graph_matrix = new int64_t*[n_vertex];
    for(int i = 0; i < n_vertex; ++i)
        graph_matrix[i] = new int64_t[n_vertex];

    for(int i = 0; i < n_vertex; ++i)
        for(int j = 0; j < n_vertex; ++j)
            graph_matrix[i][j] = 0;
}

Graph::~Graph() {
    delete[] graph_matrix;
}

void Graph::AddEdge(size_t from, size_t to) {
    try {

        if ( (from >= vertex_quan) || (to >= vertex_quan) )
            throw vertex_quan;

        graph_matrix[from][to] = rand() % 5;
//        graph_matrix[to][from] = 1;
    } catch (size_t vertex_quan) {
        std::cout << "\nInput edge ("<< from << ", " << to << ") doesn't exist in the graph of " << vertex_quan << std::endl;
        return;
    }

}

void CreateRandomGraph(Graph &g, int dense, bool PrintGraph=false) {
    int mul_size = 3 + rand() % dense;
    int edge_number = 0;
    for (int i = 0; i < g.vertex_quan * mul_size; i++) {
        int from = rand() % g.vertex_quan;
        int to = rand() % g.vertex_quan;
        if (from == to)
            continue;
        g.AddEdge(from, to);
        edge_number++;
    }
    std::cout << "Graph with " << g.vertex_quan << " vertexes and " << edge_number << " edges has been created" << std::endl;
    if (PrintGraph) {
        for (int i = 0; i < g.vertex_quan; i++) {
            for (int j = 0; j < g.vertex_quan; j++)
                std::cout << g.graph_matrix[i][j] << ' ';
            std::cout << std::endl;
        }
    }
}

// С использованием описанной выше структуры реализуйте алгоритм Дейкстры. Сгенерируйте случайные
// большие графы и оцените скорость работы на них при использовании 2, 4 и 8 - кучи.
// Отдельно рассмотрите случаи плотных и разреженных графов.

// Используйте интерфейс графа, описанный в первом бонусном задании или напишите свой.

std::vector<int64_t> DeijkstraAlgorithm(Graph &graph, int n_vertex, int d, int start_from);

int main () {

    int n_vertex = 7;
    int d = 2;
    int start_from = 0;
// Корректность__________________________________________________________________##########
    Graph graph(n_vertex);
    CreateRandomGraph(graph, n_vertex / 2, true);

    std::vector<int64_t> p = DeijkstraAlgorithm(graph, n_vertex, d, start_from);

    std::cout << std::endl;
    for(auto i = p.begin(); i != p.end(); i++)
        std::cout << *i << ' ';
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
// ______________________________________________________________________________##########


// Случайный плотный граф________________________________________________________##########

    n_vertex = 5000;

    Graph danse_graph(n_vertex);
    CreateRandomGraph(danse_graph, n_vertex);

    d = 2;
    float fTimeStart = clock()/(float)CLOCKS_PER_SEC;

    DeijkstraAlgorithm(danse_graph, n_vertex, d, start_from);

    float fTimeStop = clock()/(float)CLOCKS_PER_SEC;
    printf("\nДлительность выполнения Дейксры с %d-кучей для плотного графа -- %f секунд\n", d, fTimeStop-fTimeStart);

    d = 4;
    fTimeStart = clock()/(float)CLOCKS_PER_SEC;

    DeijkstraAlgorithm(danse_graph, n_vertex, d, start_from);

    fTimeStop = clock()/(float)CLOCKS_PER_SEC;
    printf("\nДлительность выполнения Дейксры с %d-кучей для плотного графа -- %f секунд\n", d, fTimeStop-fTimeStart);

    d = 8;
    fTimeStart = clock()/(float)CLOCKS_PER_SEC;

    DeijkstraAlgorithm(danse_graph, n_vertex, d, start_from);

    fTimeStop = clock()/(float)CLOCKS_PER_SEC;
    printf("\nДлительность выполнения Дейксры с %d-кучей для плотного графа -- %f секунд\n", d, fTimeStop-fTimeStart);
    std::cout << std::endl;
    std::cout << std::endl;


// ______________________________________________________________________________##########

// Случайный разреженный граф____________________________________________________##########

    n_vertex = 1500;

    Graph sparce_graph(n_vertex);
    CreateRandomGraph(sparce_graph, n_vertex / 10);

    d = 2;
    fTimeStart = clock()/(float)CLOCKS_PER_SEC;

    DeijkstraAlgorithm(sparce_graph, n_vertex, d, start_from);

    fTimeStop = clock()/(float)CLOCKS_PER_SEC;
    printf("\nДлительность выполнения Дейксры с %d-кучей для разреженного графа -- %f секунд\n", d, fTimeStop-fTimeStart);

    d = 4;
    fTimeStart = clock()/(float)CLOCKS_PER_SEC;

    DeijkstraAlgorithm(sparce_graph, n_vertex, d, start_from);

    fTimeStop = clock()/(float)CLOCKS_PER_SEC;
    printf("\nДлительность выполнения Дейксры с %d-кучей для разреженного графа -- %f секунд\n", d, fTimeStop-fTimeStart);

    d = 8;
    fTimeStart = clock()/(float)CLOCKS_PER_SEC;

    DeijkstraAlgorithm(sparce_graph, n_vertex, d, start_from);

    fTimeStop = clock()/(float)CLOCKS_PER_SEC;
    printf("\nДлительность выполнения Дейксры с %d-кучей для разреженного графа -- %f секунд\n", d, fTimeStop-fTimeStart);
    std::cout << std::endl;
    std::cout << std::endl;


// ______________________________________________________________________________##########

// Свои данные___________________________________________________________________##########


    std::cout << "\nВведите кол-во вершин в графе: ";
    std::cin >> n_vertex;
    std::cout << "\nВведите количество потомков в куче: ";
    std::cin >> d;
    std::cout << "\nВведите стартовую вершину: ";
    std::cin >> start_from;


    Graph my_graph(n_vertex);
    CreateRandomGraph(my_graph, n_vertex / 2, true);

    p = DeijkstraAlgorithm(my_graph, n_vertex, d, start_from);

    std::cout << std::endl;
    for(auto i = p.begin(); i != p.end(); i++)
        std::cout << *i << ' ';

    return 0;
}



std::vector<int64_t> DeijkstraAlgorithm(Graph &graph, int n_vertex, int d, int start_from) {

    std::vector<int64_t> p (n_vertex, std::numeric_limits<int64_t>::max());
    p[start_from] = 0;

    DijkstraStorage q(n_vertex, d);
    q.SetDistance(start_from, 0);
    q.MakeVertexActive(start_from);

    while(!q.Empty()) {
        int cur_d = q.GetDistanceToVertex(q.GetMinIndex()), v = q.RemoveMin();

        if (cur_d > p[v] || cur_d == -1)
            continue;

        for (size_t i = 0; i < graph.vertex_quan; i++) {
            int len = graph.graph_matrix[i][v];
            if (len != 0) {
                if (p[v]+len < p[i]) {
                    p[i] = p[v]+len;
                    q.SetDistance(i, p[i]);
                    q.MakeVertexActive(i);
                }
            }
        }
    }

    return p;
}