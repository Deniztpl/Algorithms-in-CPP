#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <limits>
#include <algorithm>
#include <fstream>

using namespace std;

// An attempt to model graphs
class Graph {
private:
    struct Edge {
        bool connected;
        float weight;
    };

    struct Vertex {
        int vertexNumber;
        float vertexValue;
    };

    int size;           // Size of graph, number of vertices
    Vertex* vertices;   // 1D array of vertices
    Edge** graph;       // 2D array of edges

public:
    // Initialize a constructor that constructs an empty graph with given size.
    Graph(int size) : size(size) {
        graph = new Edge*[size];
        vertices = new Vertex[size];
        for (int i = 0; i < size; ++i) {
            vertices[i].vertexNumber = i;
            vertices[i].vertexValue = 0;
            graph[i] = new Edge[size];
            for (int j = 0; j < size; ++j) {
                graph[i][j].connected = false;
                graph[i][j].weight = 0;
            }
        }
    }

    Graph(const std::string& filename) {
        std::ifstream file(filename);

        if (!file.is_open()) {
            std::cerr << "Error opening the file: " << filename << std::endl;
            return;
        }

        file >> size;

        // Allocate memory for the graph
        graph = new Edge*[size];
        vertices = new Vertex[size];

        for (int i = 0; i < size; ++i) {
            vertices[i].vertexNumber = i;
            vertices[i].vertexValue = 0;
            graph[i] = new Edge[size];
            for (int j = 0; j < size; ++j) {
                graph[i][j].connected = false;
                graph[i][j].weight = 0;
            }
        }

        int from, to, cost;

        while (file >> from >> to >> cost) {
            if (from >= 0 && from < size && to >= 0 && to < size) {
                graph[from][to].connected = true;
                graph[to][from].connected = true;  // Make it undirected
                graph[from][to].weight = cost;
                graph[to][from].weight = cost;  // Make it undirected
            } else {
                std::cerr << "Invalid vertex indices in the file." << std::endl;
                return;
            }
        }

        file.close();
    }



    // Deconstructor.
    ~Graph() {
        for (int i = 0; i < size; ++i)
            delete[] graph[i];
        delete[] graph;
        delete[] vertices;
    }

    // Return number of vertices.
    int V() {
        return size;
    }

    // Return number of edges.
    int E() {
        int edgeCount = 0;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (graph[i][j].connected == true) {
                    edgeCount++;
                }
            }
        }
        return (edgeCount / 2);     // Since graph is undirected
    }

    // Check if node x and node y adjacent.
    bool adjacent(int x, int y) {
        if (x >= 0 && x < V() && y >= 0 && y < V()) {
            return graph[x][y].connected;
        } else {
            return false;
        }
    }

    // List neighbour nodes of x.
    void neighbors(int x) {
        cout << "There are edges from vertex " << x << " to these vertices: " << endl;
        for (int i = 0; i < size; i++) {
            if (graph[x][i].connected == true) {
                cout << i << endl;
            }
        }
    }

    // Add edge between node i and node j.
    void addEdge(int i, int j) {
        if (i >= 0 && i < V() && j >= 0 && j < V()) {
            graph[i][j].connected = true;
            graph[j][i].connected = true;   // Make it undirected
        }
    }

    // Delete the edge between node i and node j.
    void deleteEdge(int i, int j) {
        if (i >= 0 && i < V() && j >= 0 && j < V()) {
            graph[i][j].connected = false;
            graph[j][i].connected = false;
        }
    }

    // Set and get functions for edges and vertices.
    void set_vertex_value(int x, float a) {
        vertices[x].vertexValue = a;
    }

    float get_vertex_value(int x) {
        return vertices[x].vertexValue;
    }

    void set_edge_value(int i, int j, float v) {
        if (i >= 0 && i < V() && j >= 0 && j < V()) {
            graph[i][j].weight = v;
            graph[j][i].weight = v;
        }
    }

    float get_edge_value(int i, int j) {
        return graph[i][j].weight;
    }

    // Fill a graph with random density and weight range.
    void fillWithRandom(float density, float minWeight, float maxWeight) {
        srand(time(0));
        for (int i = 0; i < size; ++i) {
            for (int j = i + 1; j < size; ++j) {
                // If random number between 0 and 1.0 less than density execute
                if ((static_cast<float>(rand()) / RAND_MAX) < density) {
                    float weight = static_cast<float>(rand()) / RAND_MAX * (maxWeight - minWeight) + minWeight;     // Random weight
                    addEdge(i, j);
                    set_edge_value(i, j, weight);
                }
            }
        }
    }

    // Print the graph as 2D array of boolean values
    void print() {
        cout << "Adjacency Matrix:" << endl;
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                cout << graph[i][j].connected << " ";
            }
            cout << endl;
        }
    }
};


// A priority queue model to be used in Dijkstra algorithm
class PriorityQueue {
private:
    struct QueueElement {
        int element;
        float priority;

        QueueElement(int el, float prio) : element(el), priority(prio) {}
    };

    vector<QueueElement> pq;

public:
    // A basic function to be used in sorting queue elements.
    static bool compareQueueElements(const QueueElement& a, const QueueElement& b) {
        return a.priority > b.priority;
    }

    // Change priority of an element.
    void chgPriority(int element, float priority) {
        for (int i = 0; i < pq.size(); ++i) {
            if (pq[i].element == element) {
                pq[i].priority = priority;
                sort(pq.begin(), pq.end(), compareQueueElements);
                return;
            }
        }
    }

    // Pop the smallest priority element from queue and return it
    int minPriority() {
        if (!pq.empty()) {
            int topElement = pq.back().element;
            pq.pop_back();
            return topElement;
        }
        return -1;
    }

    // Check if element is in queue.
    bool contains(int element) {
        for (const QueueElement& e : pq) {
            if (e.element == element) {
                return true;
            }
        }
        return false;
    }

    // Insert an element into queue according to its priority.
    void insert(int element, float priority) {
        pq.push_back(QueueElement(element, priority));
        sort(pq.begin(), pq.end(), compareQueueElements);
    }

    // Return top element(no pop).
    int top() {
        if (!pq.empty()) {
            return pq.back().element;
        }
        return -1;
    }

    // Return size of queue
    int size() {
        return pq.size();
    }
};


// A shortest path model that uses classes Graph and PriorityQueue
class ShortestPath {
private:
    Graph graph;        // class Graph
    PriorityQueue pq;   // class PriorityQueue

public:
    // Constructor
    ShortestPath(const Graph& g) : graph(g) {}

    // Find the shortest path between vertices u and w and return the sequence as vector of integers.
    vector<int> path(int u, int w) {
        int V = graph.V();
        vector<float> dist(V, numeric_limits<float>::max());    // Set all distances to vertex u to infinity
        vector<int> prev(V, -1);                                // Store previous vertices for each node.
        vector<bool> visited(V, false);                         // Store visited vertices

        dist[u] = 0;    // Self distance initialized to 0
        pq.insert(u, 0);

        // Until all nodes visited in the priority queue
        while (pq.size() > 0) {
            int current = pq.minPriority();
            visited[current] = true;

            // Check for adjacent nodes of current node. If there is a shorter way to go to node v update the priority queue.
            for (int v = 0; v < V; ++v) {
                if (graph.adjacent(current, v) && !visited[v]) {
                    float weight = dist[current] + graph.get_edge_value(current, v);
                    if (weight < dist[v]) {
                        dist[v] = weight;
                        prev[v] = current;
                        if (pq.contains(v)) {
                            pq.chgPriority(v, weight);
                        } else {
                            pq.insert(v, weight);
                        }
                    }
                }
            }
        }

        // Build the path sequence
        vector<int> pathSequence;
        int currentVertex = w;      // From node w to node u traceback
        while (currentVertex != -1) {
            pathSequence.push_back(currentVertex);
            currentVertex = prev[currentVertex];
        }

        reverse(pathSequence.begin(), pathSequence.end());  // Reverse the vector so the first element will be u
        return pathSequence;
    }

    // Return the distance between vertex u and vertex w.
    float path_size(int u, int w) {
        vector<int> pathSequence = path(u, w);
        float pathCost = 0.0;

        if (!pathSequence.empty()) {
            for (int i = 0; i < pathSequence.size() - 1; ++i) {
                int currentVertex = pathSequence[i];
                int nextVertex = pathSequence[i + 1];
                pathCost += graph.get_edge_value(currentVertex, nextVertex);
            }
        }
        return pathCost;
    }
};


// A function to print average path length.
void averageDistance (ShortestPath &shortestPath, int numVertices) {
    int startVertex = 0;
    float totalDistance = 0;      // To calculate average
    for (int i = 0; i < 50; i++){
        int finishVertex = i;
        float pathDistance = shortestPath.path_size(startVertex, finishVertex);
        totalDistance += pathDistance;

        cout << "Distance from vertex " << startVertex << " to " << finishVertex << ": " << pathDistance << endl;
    }

    float averageDistance = totalDistance / (numVertices - 1);   // Average distance of 49 possible path

    cout << endl;
    cout << "Average path length is: " << averageDistance << endl;
    cout << endl;
}


int main() {
    int numVertices = 50;
    Graph graph1(numVertices);
    graph1.fillWithRandom(0.2, 1.0, 10.0); // Fill the graph with density 0.2
    ShortestPath shortestPath1(graph1);

    Graph graph2(numVertices);  // Another instance
    graph2.fillWithRandom(0.4, 1.0, 10.0); // Fill the graph with density 0.4
    ShortestPath shortestPath2(graph2);

    cout << "Distances and average path length on graph with density 0.2 and size 50" << endl;
    cout << "-----------------------------------------------------------------------" << endl;
    averageDistance(shortestPath1, numVertices);

    cout << "******\n" << endl;

    cout << "Distances and average path length on graph with density 0.4 and size 50" << endl;
    cout << "-----------------------------------------------------------------------" << endl;
    averageDistance(shortestPath2, numVertices);

    Graph graphaa("sample_data.txt");
    graphaa.print();

     if (__cplusplus == 202002L) std::cout << "C++20\n";
    else if (__cplusplus == 201703L) std::cout << "C++17\n";
    else if (__cplusplus == 201402L) std::cout << "C++14\n";
    else if (__cplusplus == 201103L) std::cout << "C++11\n";
    else if (__cplusplus == 199711L) std::cout << "C++98\n";
    else std::cout << "pre-standard C++\n";

    return 0;
}
