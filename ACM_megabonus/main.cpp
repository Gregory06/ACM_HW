/*
Данная задача является бонусной высокого уровня, она гарантирует вам +1.5 балла при чистом решении к 10-бальной оценке в конце
семестра в случае, если ваша оценка на тот момент будет >= уд(3).

Тем не менее дедлайн для отправки решений - 16 марта, 23:00.

Решать данную задачу на интерпретируемых языках не имеет смысла, так как они дадут просадку производительности "из коробки"
Желательно - С++, но готов рассматривать Go, Rust, C.

Неполные решения (например реализация лишь одного из алгоритмов, будут так же оценены, возможно в +0.5 -- +1 балл к итоговой)
*/

/*
    Хотя в нашем курсе класс задач и языков NP рассматривается как нечто трудное,
    что не стоит даже пытаться решать за полином, на практике подобные задачи решать все же приходится.
    Обычно используются различные эвристики, которые все так же плохи в худшем случае, как и полный перебор,
    но способные на большей части возможных входов давать быстрый ответ.

    В этот раз мы предлагаем вам решить задачу поиска клик в неориентированных графах без петель и кратных ребер,
    используя 2 различных алгоритма, а потом сравнить их результат.

    1. Детерминированный алгоритм, базирующийся на принципе meet-in-the-middle. На семинаре вам будет рассказано
        применение этой идеи для решения NP-полной задачи, несколько похожей на SUBSET-SUM (видимо, пока не будет, поэтому конспект в конце файла),
        а так же задачи дискретного логарифмирования: https://en.wikipedia.org/wiki/Baby-step_giant-step, он же -- Алгоритм Шенкса

        Вам стоит самостоятельно реализовать проверку работы реализации и замер времени работы. От вас ожидается, что граф на
        48 вершинах будет обрабатываться примерно секунду. (Мне с помощью алгоритмических и технических трюков, основанных
        на особенностях работы процессоров intel удалось обрабатывать 54 вершины за 0.8 секунд, сложность имеет вид ~2^(N/2) * N)
        Этот алгоритм показывает примерно одинаковое время работы на любых графах и не имеет сильно выраженных лучших и худших случаев.

    2. Эвристический алгоритм "Ветвей и границ", он же "Branches and bounds". Вот полезная статья на эту тему: https://bit.ly/2VwVGta.
        Эвристические алгоритмы позволяют осуществлять полный перебор с выбором некоторых приоритетных направлений, которые
        скорее всего и будут содержать решение. Не без ложки дегтя: для большинства эвристик существуют контрпримеры, для
        которых они показывают просто ужасные результаты, хоть и решают правильно.

        В практических задачах обычно сначала используется этот подход. В случае, если он не дает решения за некоторое время
        происходит смена алгоритма на другой.
*/

/*
    Граф стоит хранить в отдельном классе для обеспечения архитектуры и читаемости кода
*/

#include <vector>
#include <cstdint>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <time.h>
#include <cstdlib>

bool comparePairs(std::pair<unsigned long long ,bool> a, std::pair<unsigned long long ,bool> b);
bool comparePairAndInt(std::pair<unsigned long long ,bool> a, unsigned long long b);

class Graph {
public:
    int **graph_matrix;
    size_t vertex_quan;
    Graph(size_t n_vertex);
    ~Graph();
    void AddEdge(size_t from, size_t to);

    // Метод возвращает true, если в графе действительно есть клика размера заданного размера
    virtual bool HasClique(size_t clique_size) = 0;

    // Метод возвращает размер максимальной клики для графа
    virtual size_t GetMaxCliqueSize() = 0;

    // Метод возвращает вершины максимальной клики
    virtual std::vector<size_t> GetMaxClique() = 0;
};

Graph::Graph(size_t n_vertex) {
    vertex_quan = n_vertex;
    graph_matrix = new int*[n_vertex];
    for(int i = 0; i < n_vertex; ++i)
        graph_matrix[i] = new int[n_vertex];

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

        graph_matrix[from][to] = 1;
        graph_matrix[to][from] = 1;
    } catch (size_t vertex_quan) {
        std::cout << "\nInput edge ("<< from << ", " << to << ") doesn't exist in the graph of " << vertex_quan << std::endl;
        return;
    }

}

/*
    От приведенного выше класса вам стоит унаследовать 2 потомков, первый будет предоставлять реализацию MITM
    алгоритма, а второй - эвристики.

    Определить способ хранения можно как в обоих потомках (возможно различным), но реализацию AddEdge и хранения графа
    я бы предложит реализовывать в Graph.

    На любые вопросы по С++ готов ответить в личке или после семинаров.
*/

struct Comp {
    bool operator() (const std::pair<unsigned long long, bool> a, unsigned long long b);

    bool operator() (unsigned long long b, const std::pair<unsigned long long, bool> a);
};

class GraphMITM : public Graph {
public:
    GraphMITM(size_t n_vertex): Graph(n_vertex) {}

    bool CheckPrevVertexes(size_t new_vertex, int other_vertex) {
        for (int i = 0; i < new_vertex; i++) {
            if (bool((1 << i)  &  other_vertex))
                if (!graph_matrix[new_vertex][i])
                    return false;
        }
        return true;
    }

    std::vector<std::pair<unsigned long long ,bool>> SubGraphProcess(size_t vertex_from, size_t vertex_to) {
        std::vector<std::pair<unsigned long long ,bool>> is_clique;
        is_clique.push_back(std::pair(0 ,true));
        for(unsigned long long i = 1; i < pow(2, vertex_to - vertex_from); i++) {
            unsigned long long current = i * pow(2, vertex_from);
            int prev = i;
            size_t new_vertex = size_t (floor(log2(float(i * pow(2, vertex_from)))));
            current &= ~(1 << new_vertex);
            prev &= ~(1 << size_t (floor(log2(float(i)))));
            bool is_in_clique = CheckPrevVertexes(new_vertex, current);
            if ( is_clique[prev].second && is_in_clique )
                is_clique.push_back(std::pair((unsigned long long)(i * pow(2, vertex_from)) ,true));
            else
                is_clique.push_back(std::pair((unsigned long long)(i * pow(2, vertex_from)) ,false));
        }
        return is_clique;
    }

    bool IsOneClique(std::pair<unsigned long long, bool> a, std::pair<unsigned long long, bool> b, int half) {
        int logf = int (floor(log2(float(a.first)))) + 1;
        int logs = int (floor(log2(float(b.first)))) + 1;
        unsigned long long c = a.first * pow(2, 0);

        for (int i = 0; i < logf; i++)
            if (bool((1 << i)  &  a.first))
                for (int j = half; j < logs; j++) {
                    if (bool((1 << j) &  b.first)) {
                        if (!graph_matrix[i][j])
                            return false;
                    }
            }
        return true;
    }

    bool HasCliqueSlow(size_t clique_size) {
        std::vector<std::pair<unsigned long long, bool>> combinations;

        combinations = SubGraphProcess(0, vertex_quan);
        std::sort(combinations.begin(), combinations.end(), comparePairs);

        for (int i = 0; i < combinations.size(); ++i) {
            if (combinations[i].second) {
                int vert_num = 0, a = combinations[i].first;
                while (a > 0) {
                    if (a % 2 != 0) {
                        vert_num++;
                    }
                    a /= 2;
                }
                int k = clique_size - vert_num;

                if (k == 0) {
                    unsigned long long c = combinations[i].first * pow(2, 0);
                    for (int j = sizeof(unsigned long long) * 8 - 1; j > -1; --j) {
                        std::cout << ((c >> j) & 1);
                    }

                    return true;
                }
            }
        }
        std::cout << "\nNONE" << std::endl;
        return false;
    }


    virtual bool HasClique(size_t clique_size) {
        int half = int(vertex_quan / 2);
        std::vector<std::pair<unsigned long long ,bool>> first_half, second_half;

        first_half = SubGraphProcess(0, half);
        second_half = SubGraphProcess(half, vertex_quan);

        std::sort(second_half.begin(), second_half.end(), comparePairs);
        for (int i = 0; i < first_half.size(); ++i) {
            if (first_half[i].second) {
                int vert_num = 0;
                unsigned long long a = first_half[i].first;
                while (a > 0)
                {
                    if (a % 2 != 0)
                    {
                        vert_num++;
                    }
                    a /= 2;
                }
                unsigned long long k = clique_size - vert_num;

                if (k == 0) {

                    unsigned long long c = first_half[i].first * pow(2, 0);
                    for (int j = sizeof(unsigned long long) * 8 - 1; j > -1; --j) {
                        std::cout << ((c >> j) & 1);
                    }

                    return true;

                } else {
                    auto pairs =  std::equal_range(second_half.begin(), second_half.end(), k, Comp{});
                    for (auto pair = pairs.first; pair != pairs.second; pair++) {
                        if (pair->second) {
                            if (IsOneClique(first_half[i], *pair, half)) {
                                unsigned long long c = first_half[i].first + pair->first;
                                for (int j = sizeof(unsigned long long) * 8 - 1; j > -1; --j) {
                                    std::cout << ((c >> j) & 1);
                                }
                                return true;
                            }
                        }
                    }
                }
            }
        }
        std::cout << "\nNONE" << std::endl;
        return false;
    }

    virtual size_t GetMaxCliqueSize() {

    }

    virtual std::vector<size_t> GetMaxClique() {

    }

};

bool comparePairs(std::pair<unsigned long long ,bool> a, std::pair<unsigned long long ,bool> b) {
    unsigned long long f = a.first;
    unsigned long long s = b.first;
    int countf = 0, counts = 0;
    while (f > 0)
    {
        if (f % 2 != 0)
        {
            countf++;
        }
        f /= 2;
    }
    while (s > 0)
    {
        if (s % 2 != 0)
        {
            counts++;
        }
        s /= 2;
    }
    return (countf < counts);
}

bool Comp::operator() (const std::pair<unsigned long long ,bool> a, unsigned long long b) {
    unsigned long long f = a.first;
    int countf = 0;
    while (f > 0)
    {
        if (f % 2 != 0)
        {
            countf++;
        }
        f /= 2;
    }
    return (countf < b);
}

bool Comp::operator() ( unsigned long long b, const std::pair<unsigned long long ,bool> a) {
    unsigned long long f = a.first;
    int countf = 0;
    while (f > 0)
    {
        if (f % 2 != 0)
        {
            countf++;
        }
        f /= 2;
    }
    return (b < countf);
}

class GraphBRAB : public Graph {
    // ..
};


/*
    Задача, демонстрирующая еще один способ применения идеи MITM

    Условие:
        На вход программе подаются N целых чисел, нужно найти два непересекающихся непустых подмножества этих чисел таких,
        что суммы чисел каждого подмножества будут равны.

    Решение в лоб:
        Рассмотрим все N-значные числа в троичной системе исчисления, в каждом из них 0 на i позиции будет означать, что i число не взять,
        1 - что оно принадлежит первому подмножеству, 2 - второму. Для каждого числа проверим условие, что соответствующие подмножества не пусты и
        их суммы совпадают. Оценка сложности снизу -- \Omega(3^N), как перебор количества разбиений. На самом деле все еще немного хуже,
        так как для каждого разбиения нужно считать сумму.

    Решение по MITM:
        1. Разобьем числа на две группы по N/2. (если N нечетно, то группы будут отличаться размерном на 1)
        2. Для каждой группы построим аналогично предыдущему решению 3^(N/2) разбиений и вычислим их суммы, причем
            в этот раз включим в перебор и такие разбиения, где одно из множеств (а может и оба) будет пустым.
            Сохраним не обе суммы, а лишь только разность между ними вместе с описанием разбиений в массивы A и B.
            (А хранит пары <разбиение, разность сумм подмножеств> для первой половины всех чисел, B - для второй,
            итого их размеры в районе 3^(N/2))

        3. Отсортируем A и B по возрастанию и за линейное по их сумме длин время будем перебирать все пары a[i] \in A, b[i] \in B, таких,
            что a[i] == b[i].

        4. Чуть-чуть подумаем над процессом из п. 3 и выведем ответ.

        Идея: Если есть два равных числа в массивах, то a[l] = \sum_i A[s_i] - \sum_j A[k_j], то b[q] = \sum_r B[p_r] - \sum_t B[u_t],
            Тогда \sum_i A[s_i] + \sum_t B[u_t] = \sum_r B[p_r] + \sum_j A[k_j] - мы нашли два равных по сумме непересекающихся подмножества.
            Проверка на непустоту осуществляется отдельно.

*/

void CreateRandomGraph(Graph &g);
void CreateBadGraph(Graph& result);

int main() {

    size_t n_vertex = 0, clique_size = 0;

    std::cout << "\nВведите кол-во вершин в графе: ";
    std::cin >> n_vertex;

    std::cout << "\nВведите кол-во вершин в клике: ";
    std::cin >> clique_size;

    GraphMITM graph(n_vertex);

    CreateBadGraph(graph);

    for (int i = 0; i < graph.vertex_quan; i++) {
        for (int j = 0; j < graph.vertex_quan; j++)
            std::cout << graph.graph_matrix[i][j] << ' ';
        std::cout << std::endl;
    }


    float fTimeStart = clock()/(float)CLOCKS_PER_SEC;

    graph.HasClique(clique_size);

    float fTimeStop = clock()/(float)CLOCKS_PER_SEC;
    printf("\nДлительность выполнения meet-in-the-middle %f секунд\n", fTimeStop-fTimeStart);

    fTimeStart = clock()/(float)CLOCKS_PER_SEC;

    graph.HasCliqueSlow(clique_size);

    fTimeStop = clock()/(float)CLOCKS_PER_SEC;
    printf("\nДлительность выполнения стандартного перебора %f секунд\n", fTimeStop-fTimeStart);
}

void CreateRandomGraph(Graph &g) {
    int mul_size = 3 + rand() % g.vertex_quan;
    for (int i = 0; i < g.vertex_quan * mul_size; i++) {
        int from = rand() % g.vertex_quan;
        int to = rand() % g.vertex_quan;
        g.AddEdge(from, to);
    }
}

void CreateBadGraph(Graph& result) {
    size_t n_vertices = result.vertex_quan;
    for (int i = 0; i < n_vertices; ++i) {
        for (int j = 0; j < n_vertices; ++j) {
            if (i >= j) {
                continue;
            }
            if ((i < 5) && (j > 35)) {
                continue;
            } else {
                result.AddEdge(i, j);
            }
        }
    }
}