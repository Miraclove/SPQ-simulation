#include <iostream>
#include <queue>
#include <vector>
#include <fstream>
#include <sstream>
#include <random>
#include <cmath>
#include <cstdlib>
#include <iomanip>

using namespace std;

// Constants
const double SERVICE_RATE = 10e6; // 10 Mbps (in bits/sec)

// Enumeration for traffic type / priority
enum TrafficType { AUDIO, VIDEO, DATA };

// A packet structure
struct Packet {
    double arrivalTime;  // time when the packet arrived at the node
    double sizeBits;     // packet size in bits
    TrafficType type;    // determines its priority (AUDIO highest, then VIDEO, then DATA)
    bool isReference;    // flag if this packet belongs to the reference flow
    int hop;             // current node number in the path (1...M)
    double origTime;     // time when the reference packet was generated (for end-to-end delay)
};

// Event types
enum EventType { ARRIVAL, DEPARTURE };

// An event structure for the simulation
struct Event {
    double time;       // simulation time at which the event occurs
    EventType type;    // ARRIVAL or DEPARTURE
    int nodeID;        // node where the event occurs (1...M)
    Packet packet;     // the packet associated with the event

    // For priority_queue (min-heap on time)
    bool operator<(const Event &other) const {
        return time > other.time; // reverse since priority_queue is max-heap by default
    }
};

// Structure to hold node statistics (per priority queue)
struct NodeStats {
    // For each priority (index 0: AUDIO, 1: VIDEO, 2: DATA)
    double totalDelay[3] = {0, 0, 0};
    int servedCount[3] = {0, 0, 0};
    int dropCount[3] = {0, 0, 0};
    int arrivalCount[3] = {0, 0, 0};

    // For time-average queue length (area under queue curve)
    double areaQueue[3] = {0, 0, 0};
    // Last update time for each node
    double lastUpdateTime = 0;
};

// A node structure. Each node has three FIFO queues (one per priority) and a busy flag.
struct Node {
    int id;
    // queues for each priority: AUDIO (index 0), VIDEO (index 1), DATA (index 2)
    queue<Packet> queues[3];
    bool busy = false;
    int capacity;  // total capacity (K) for all packets in the node (waiting + in service); negative means infinite
    NodeStats stats;
};

// Global simulation parameters (set from command line)
int M;                      // number of nodes
int K;                      // queue capacity per node (total for three queues; negative means infinite)
double simDuration;         // simulation duration in seconds
string refTypeStr;          // reference type ("audio", "video", "data")
double loadMin, loadMax, loadStep; // offered load range (ρ values)
bool debugMode = false;     // Flag for debug output

// Traffic source parameters structure
struct TrafficParams {
    double peakKbps;      // peak rate (kbps)
    double meanOn;        // average ON time (sec)
    double meanOff;       // average OFF time (sec)
    int packetSizeBytes;  // nominal packet size in bytes
};

// Fixed parameters from the assignment table:
const TrafficParams AUDIO_PARAMS = {64, 0.36, 0.64, 120};
const TrafficParams VIDEO_PARAMS = {384, 0.33, 0.73, 1000};
const TrafficParams DATA_PARAMS  = {256, 0.35, 0.65, 583};

// Pre-compute average packet generation rates per source (packets/sec)
double avgRatePerSource(TrafficParams p) {
    double onFraction = p.meanOn / (p.meanOn + p.meanOff);
    double packetTime = (p.packetSizeBytes * 8) / (p.peakKbps * 1000);
    return onFraction / packetTime;
}

double AUDIO_RATE = avgRatePerSource(AUDIO_PARAMS);
double VIDEO_RATE = avgRatePerSource(VIDEO_PARAMS);
double DATA_RATE  = avgRatePerSource(DATA_PARAMS);

// Average offered load (in bits/sec) per source
double offeredLoadPerSource(TrafficParams p) {
    double onFraction = p.meanOn / (p.meanOn + p.meanOff);
    return p.peakKbps * 1000 * onFraction; // Convert kbps to bps and multiply by ON fraction
}

double AUDIO_LOAD = offeredLoadPerSource(AUDIO_PARAMS);
double VIDEO_LOAD = offeredLoadPerSource(VIDEO_PARAMS);
double DATA_LOAD  = offeredLoadPerSource(DATA_PARAMS);

// Random number generator
default_random_engine generator((unsigned)time(nullptr));

// Exponential random variate generator
double expRandom(double lambda) {
    exponential_distribution<double> dist(lambda);
    return dist(generator);
}

// Function to schedule the next arrival event for a given node and traffic type
void scheduleArrivalEvent(priority_queue<Event> &eventQueue, int nodeID, TrafficType t, double currentTime, double arrivalRate, const TrafficParams &params, bool isRef = false) {
    if(arrivalRate <= 0) return;
    double dt = expRandom(arrivalRate);
    double eventTime = currentTime + dt;
    Packet p;
    p.arrivalTime = eventTime;
    p.sizeBits = params.packetSizeBytes * 8;
    p.type = t;
    p.isReference = isRef;
    p.hop = nodeID; 
    p.origTime = eventTime; // mark generation time for reference packets

    Event ev;
    ev.time = eventTime;
    ev.type = ARRIVAL;
    ev.nodeID = nodeID;
    ev.packet = p;
    eventQueue.push(ev);
}

// When a packet finishes service at a node, schedule its departure event.
void scheduleDepartureEvent(priority_queue<Event> &eventQueue, int nodeID, const Packet &p, double currentTime) {
    double serviceTime = p.sizeBits / SERVICE_RATE;
    Event ev;
    ev.time = currentTime + serviceTime;
    ev.type = DEPARTURE;
    ev.nodeID = nodeID;
    ev.packet = p;
    eventQueue.push(ev);
}

// Helper function: returns total number of packets in the system at a node (waiting plus in service)
int totalSystemPackets(const Node &node) {
    return node.queues[0].size() + node.queues[1].size() + node.queues[2].size() + (node.busy ? 1 : 0);
}

// Main simulation function for one offered load value.
// Returns a string (CSV-formatted lines) with the measured statistics.
string simulateOneLoad(double target_rho, TrafficType refTrafficType) {
    // Determine reference traffic parameters based on type.
    TrafficParams refParams;
    double refRate;
    if(refTrafficType == AUDIO) {
        refParams = AUDIO_PARAMS;
        refRate = AUDIO_RATE;
    } else if(refTrafficType == VIDEO) {
        refParams = VIDEO_PARAMS;
        refRate = VIDEO_RATE;
    } else {
        refParams = DATA_PARAMS;
        refRate = DATA_RATE;
    }
    double refLoad = offeredLoadPerSource(refParams);

    // FIXED: Calculate required number of sources based on the target load
    // Use the example from the PDF as a reference point
    // For example, Na=4, Nv=5, Nd=5 gives ρ=0.105 according to the PDF
    
    // Base ratios from project requirements (4:5:5)
    int baseAudio = 4, baseVideo = 5, baseData = 5;
    double baseRho = 0.105; // As per the example in the PDF
    
    // More precise calculation of number of sources needed for each target load
    // This ensures a more accurate and smoother progression of offered load
    double scaleFactor = target_rho / baseRho;
    
    // FIXED: Calculate exact source counts to achieve the target load without arbitrary jumps
    int N_audio = max(1, (int)round(baseAudio * scaleFactor));
    int N_video = max(1, (int)round(baseVideo * scaleFactor));
    int N_data = max(1, (int)round(baseData * scaleFactor));
    
    // For loads very close to 1.0, ensure we have enough traffic to create congestion
    // but avoid the arbitrary jumps that were causing the discontinuities
    // if (target_rho > 0.9) {
    //     double adjustFactor = 1.0 + (target_rho - 0.9) * 2.0; // Smooth adjustment
    //     N_audio = (int)(N_audio * adjustFactor);
    //     N_video = (int)(N_video * adjustFactor);
    //     N_data = (int)(N_data * adjustFactor);
    // }
    
    // Calculate actual offered load for reporting
    double totalLoad = refLoad + (N_audio * AUDIO_LOAD) + (N_video * VIDEO_LOAD) + (N_data * DATA_LOAD);
    double actualRho = totalLoad / SERVICE_RATE;
    
    if (debugMode) {
        cout << "Target ρ: " << fixed << setprecision(3) << target_rho 
             << ", Actual ρ: " << fixed << setprecision(3) << actualRho
             << " (N_audio=" << N_audio << ", N_video=" << N_video << ", N_data=" << N_data << ")" << endl;
    }

    // For aggregated Poisson arrival processes at each node.
    double bgAudioArrivalRate = N_audio * AUDIO_RATE;
    double bgVideoArrivalRate = N_video * VIDEO_RATE;
    double bgDataArrivalRate  = N_data * DATA_RATE;
    // Reference arrival rate (only at node 1)
    double refArrivalRate = refRate;

    // Create M nodes (using indices 1..M).
    vector<Node> nodes(M+1);
    for (int i = 1; i <= M; i++) {
        nodes[i].id = i;
        nodes[i].capacity = K; // if K is negative, treat as infinite capacity
        nodes[i].busy = false;
        nodes[i].stats.lastUpdateTime = 0;
    }

    // Statistics for reference flow (end-to-end).
    int refGenerated = 0;
    int refDropped = 0;
    double sumRefDelay = 0;
    int refServed = 0;

    // Event queue.
    priority_queue<Event> eventQueue;
    double currentTime = 0;

    // Schedule initial background arrival events for each node and type.
    for (int i = 1; i <= M; i++) {
        scheduleArrivalEvent(eventQueue, i, AUDIO, currentTime, bgAudioArrivalRate, AUDIO_PARAMS, false);
        scheduleArrivalEvent(eventQueue, i, VIDEO, currentTime, bgVideoArrivalRate, VIDEO_PARAMS, false);
        scheduleArrivalEvent(eventQueue, i, DATA, currentTime, bgDataArrivalRate, DATA_PARAMS, false);
    }
    // Schedule initial reference arrival event at node 1.
    scheduleArrivalEvent(eventQueue, 1, refTrafficType, currentTime, refArrivalRate, refParams, true);
    refGenerated++;

    // FIXED: For high loads, use longer warm-up period
    double warmupPeriod = simDuration * 0.1;
    if (target_rho > 0.8) {
        warmupPeriod = simDuration * 0.2; // Longer warm-up for high loads
    }
    
    // Variables to track system state during simulation
    double lastStatsResetTime = 0;
    double maxQueueSize = 0;
    
    // Run the simulation until simDuration.
    while (!eventQueue.empty() && currentTime < simDuration) {
        Event ev = eventQueue.top();
        eventQueue.pop();
        double dt = ev.time - currentTime;
        currentTime = ev.time;

        // Update the area under the queue length curve (for time-average backlog).
        Node &node = nodes[ev.nodeID];
        for (int i = 0; i < 3; i++) {
            int qlen = node.queues[i].size();
            node.stats.areaQueue[i] += qlen * dt;
            
            // Track maximum queue size for debugging
            if (qlen > maxQueueSize) {
                maxQueueSize = qlen;
            }
        }
        node.stats.lastUpdateTime = currentTime;

        // If we've passed the warm-up period and haven't reset stats yet
        if (currentTime > warmupPeriod && lastStatsResetTime == 0) {
            // Reset statistics for all nodes to collect only steady-state metrics
            for (int i = 1; i <= M; i++) {
                for (int j = 0; j < 3; j++) {
                    nodes[i].stats.totalDelay[j] = 0;
                    nodes[i].stats.servedCount[j] = 0;
                    nodes[i].stats.dropCount[j] = 0;
                    nodes[i].stats.arrivalCount[j] = 0;
                    nodes[i].stats.areaQueue[j] = 0;
                }
                nodes[i].stats.lastUpdateTime = currentTime;
            }
            
            // Reset reference traffic statistics
            refGenerated = 0;
            refDropped = 0;
            sumRefDelay = 0;
            refServed = 0;
            
            lastStatsResetTime = currentTime;
            
            if (debugMode) {
                cout << "Warmup complete at time " << currentTime << ". Statistics reset." << endl;
            }
        }

        if (ev.type == ARRIVAL) {
            // Arrival event at node ev.nodeID.
            int prio = (int)ev.packet.type; // 0 for AUDIO, 1 for VIDEO, 2 for DATA
            node.stats.arrivalCount[prio]++;

            // Check if the system is full (waiting + in service).
            int curTotal = totalSystemPackets(node);
            bool full = (node.capacity >= 0 && curTotal >= node.capacity);

            if (full) {
                // Drop packet and update drop statistics.
                node.stats.dropCount[prio]++;
                if(ev.packet.isReference) {
                    refDropped++;
                }
                
                if (debugMode && currentTime > simDuration * 0.5 && currentTime - lastStatsResetTime > 10) {
                    cout << "DROP at node " << ev.nodeID << ", time " << fixed << setprecision(2) << currentTime 
                         << ", prio " << prio << ", total packets: " << curTotal << "/" << node.capacity << endl;
                }
            } else {
                // If the server is idle, start service immediately.
                if (!node.busy) {
                    node.busy = true;
                    scheduleDepartureEvent(eventQueue, ev.nodeID, ev.packet, currentTime);
                } else {
                    // Otherwise, add the packet to the appropriate queue.
                    node.queues[prio].push(ev.packet);
                }
            }

            // Schedule the next arrival for the same type and same node.
            if (!ev.packet.isReference) {
                if (ev.packet.type == AUDIO)
                    scheduleArrivalEvent(eventQueue, ev.nodeID, AUDIO, currentTime, bgAudioArrivalRate, AUDIO_PARAMS, false);
                else if (ev.packet.type == VIDEO)
                    scheduleArrivalEvent(eventQueue, ev.nodeID, VIDEO, currentTime, bgVideoArrivalRate, VIDEO_PARAMS, false);
                else if (ev.packet.type == DATA)
                    scheduleArrivalEvent(eventQueue, ev.nodeID, DATA, currentTime, bgDataArrivalRate, DATA_PARAMS, false);
            } else {
                // For reference arrivals, schedule only at node 1.
                if (ev.nodeID == 1) {
                    scheduleArrivalEvent(eventQueue, 1, refTrafficType, currentTime, refArrivalRate, refParams, true);
                    refGenerated++;
                }
            }
        }
        else if (ev.type == DEPARTURE) {
            // Departure event at node ev.nodeID: the packet finished service.
            int prio = (int)ev.packet.type;
            double delay = currentTime - ev.packet.arrivalTime;
            nodes[ev.nodeID].stats.totalDelay[prio] += delay;
            nodes[ev.nodeID].stats.servedCount[prio]++;

            // For reference traffic: if not at the final node, forward the packet.
            if (ev.packet.isReference) {
                if (ev.nodeID < M) {
                    Packet newPkt = ev.packet;
                    newPkt.arrivalTime = currentTime; // immediate arrival at next node
                    newPkt.hop = ev.nodeID + 1;
                    Event newEv;
                    newEv.time = currentTime;
                    newEv.type = ARRIVAL;
                    newEv.nodeID = ev.nodeID + 1;
                    newEv.packet = newPkt;
                    eventQueue.push(newEv);
                } else {
                    // Packet reached the destination.
                    sumRefDelay += (currentTime - ev.packet.origTime);
                    refServed++;
                }
            }
            
            // Check if any packet is waiting in the queues (in priority order).
            bool found = false;
            for (int i = 0; i < 3; i++) {
                if (!node.queues[i].empty()) {
                    Packet nextPkt = node.queues[i].front();
                    node.queues[i].pop();
                    scheduleDepartureEvent(eventQueue, ev.nodeID, nextPkt, currentTime);
                    found = true;
                    break;
                }
            }
            if (!found) {
                node.busy = false;
            }
        }
        
        // Periodically output debugging information
        if (debugMode && currentTime > lastStatsResetTime && (int)(currentTime) % 100 == 0) {
            // Find the maximum queue length across all nodes and priorities
            int maxQLen = 0;
            for (int i = 1; i <= M; i++) {
                for (int p = 0; p < 3; p++) {
                    if (nodes[i].queues[p].size() > maxQLen) {
                        maxQLen = nodes[i].queues[p].size();
                    }
                }
            }
            
            if (maxQLen > 0) {
                cout << "Time " << fixed << setprecision(1) << currentTime 
                     << ", Max queue length: " << maxQLen << endl;
            }
        }
    }
    
    // Check if simulation ended early
    if (currentTime < simDuration) {
        cout << "Warning: Simulation ended at time " << currentTime 
             << " before reaching the target duration of " << simDuration << " seconds." << endl;
    }

    // Calculate effective simulation duration (excluding warm-up)
    double effectiveDuration = currentTime - lastStatsResetTime;
    if (effectiveDuration <= 0) effectiveDuration = simDuration; // Fallback if no warmup reset

    // Collect statistics and output CSV lines.
    stringstream ss;
    for (int i = 1; i <= M; i++) {
        Node &node = nodes[i];
        for (int p = 0; p < 3; p++) {
            double avgDelay = (node.stats.servedCount[p] > 0) ? (node.stats.totalDelay[p] / node.stats.servedCount[p]) : 0;
            double blockingRatio = (node.stats.arrivalCount[p] > 0) ? ((double)node.stats.dropCount[p] / node.stats.arrivalCount[p]) : 0;
            double avgBacklog = node.stats.areaQueue[p] / effectiveDuration;
            
            // Debug output for checking if packets are being dropped
            if (debugMode && (node.stats.dropCount[p] > 0 || node.stats.servedCount[p] > 100) && i == 1) {
                cout << "Node " << i << " Priority " << p << ": "
                     << "Arrived=" << node.stats.arrivalCount[p] 
                     << ", Served=" << node.stats.servedCount[p]
                     << ", Dropped=" << node.stats.dropCount[p]
                     << ", AvgDelay=" << avgDelay
                     << ", BlockingRatio=" << blockingRatio
                     << ", AvgBacklog=" << avgBacklog << endl;
            }
            
            // CSV: OfferedLoad, Node, QueueType, AvgDelay, BlockingRatio, AvgBacklog, N_audio, N_video, N_data
            ss << target_rho << "," << i << "," 
               << (p==0?"Premium":(p==1?"Assured":"BestEffort")) << ","
               << avgDelay << "," << blockingRatio << "," << avgBacklog << ","
               << N_audio << "," << N_video << "," << N_data << "\n";
        }
    }
    // Output overall reference traffic statistics.
    double refDropRatio = (refGenerated > 0) ? ((double)refDropped / refGenerated) : 0;
    double avgEndToEndDelay = (refServed > 0) ? (sumRefDelay / refServed) : 0;
    
    if (debugMode) {
        cout << "Reference traffic: Generated=" << refGenerated 
             << ", Dropped=" << refDropped 
             << ", Served=" << refServed
             << ", AvgE2EDelay=" << avgEndToEndDelay
             << ", BlockingRatio=" << refDropRatio << endl;
    }
    
    ss << target_rho << ",Reference,EndToEndDelay," << avgEndToEndDelay << "\n";
    ss << target_rho << ",Reference,BlockingRatio," << refDropRatio << "\n";

    return ss.str();
}

int main(int argc, char *argv[]) {
    if(argc < 7) {
        cerr << "Usage: " << argv[0] << " M K simDuration refType loadMin loadMax loadStep(optional)" << endl;
        cerr << "Example: " << argv[0] << " 5 100 100 audio 0.1 0.9 0.1" << endl;
        return 1;
    }
    M = stoi(argv[1]);
    K = stoi(argv[2]); // use negative for infinite capacity
    simDuration = stod(argv[3]);
    refTypeStr = argv[4];
    loadMin = stod(argv[5]);
    loadMax = stod(argv[6]);
    if(argc >= 8)
        loadStep = stod(argv[7]);
    else
        loadStep = 0.1;

    // Optional argument: output file name.
    string outputFileName = "";
    if(argc >= 9) {
        outputFileName = argv[8];
    }
    
    // Check for debug flag
    if(argc >= 10 && string(argv[9]) == "debug") {
        debugMode = true;
        cout << "Debug mode enabled" << endl;
    }

    // Determine reference traffic type.
    TrafficType refTrafficType;
    if(refTypeStr == "audio")
        refTrafficType = AUDIO;
    else if(refTypeStr == "video")
        refTrafficType = VIDEO;
    else if(refTypeStr == "data")
        refTrafficType = DATA;
    else {
        cerr << "Unknown reference type. Use 'audio', 'video', or 'data'." << endl;
        return 1;
    }

    // FIXED: Increase simulation duration for high traffic loads
    if (loadMax > 0.8) {
        cout << "High traffic load detected (up to ρ=" << loadMax << "). ";
        double origDuration = simDuration;
        simDuration *= 1.5; // Increase simulation time by 50% for high loads
        cout << "Increasing simulation duration from " << origDuration << " to " << simDuration << " seconds." << endl;
    }

    // Open CSV file for writing simulation results.
    ofstream csvFile(outputFileName.empty() ? "results.csv" : outputFileName);
    if(!csvFile.is_open()) {
        cerr << "Error opening " << (outputFileName.empty() ? "results.csv" : outputFileName) << " for writing." << endl;
        return 1;
    }

    // Print simulation parameters
    cout << "Simulation Parameters:" << endl;
    cout << "Number of nodes (M): " << M << endl;
    cout << "Queue capacity (K): " << (K < 0 ? "Infinite" : to_string(K)) << endl;
    cout << "Simulation duration: " << simDuration << " seconds" << endl;
    cout << "Reference traffic type: " << refTypeStr << endl;
    cout << "Offered load range: " << loadMin << " to " << loadMax << " in steps of " << loadStep << endl;
    cout << "-----------------------------------------------" << endl;

    // Print theoretical values for reference
    cout << "Per-source average loads (bits/sec):" << endl;
    cout << "Audio: " << AUDIO_LOAD << endl;
    cout << "Video: " << VIDEO_LOAD << endl;
    cout << "Data: " << DATA_LOAD << endl;
    cout << "-----------------------------------------------" << endl;

    // Write header line with updated metric names.
    csvFile << "OfferedLoad,Node,QueueType,AvgDelay,BlockingRatio,AvgBacklog,N_audio,N_video,N_data\n";

    // Loop over offered load values.
    for(double rho = loadMin; rho <= loadMax + 1e-6; rho += loadStep) {
        cout << "Simulating offered load ρ = " << fixed << setprecision(2) << rho << "..." << endl;
        string resultLines = simulateOneLoad(rho, refTrafficType);
        csvFile << resultLines;
        cout << "Completed." << endl;
    }

    csvFile.close();
    cout << "Simulation complete. Results saved in " << (outputFileName.empty() ? "results.csv" : outputFileName) << endl;
    return 0;
}