#include <vector>
#include "../Math/Rect.h"

template<typename _Tp> static inline
double jaccardDistance(const Rectf& a, const Rectf& b) {
    _Tp Aa = a.Area();
    _Tp Ab = b.Area();

    if ((Aa + Ab) <= std::numeric_limits<_Tp>::epsilon()) {
        // jaccard_index = 1 -> distance = 0
        return 0.0;
    }

    double Aab = (a & b).Area();
    // distance = 1 - jaccard_index
    return 1.0 - Aab / (Aa + Ab - Aab);
}

template <typename T>
static inline bool SortScorePairDescend(const std::pair<float, T>& pair1,
    const std::pair<float, T>& pair2)
{
    return pair1.first > pair2.first;
}

inline void GetMaxScoreIndex(const std::vector<float>& scores, const float threshold, const int top_k,
    std::vector<std::pair<float, int> >& score_index_vec)
{
    // Generate index score pairs.
    for (size_t i = 0; i < scores.size(); ++i)
    {
        if (scores[i] > threshold)
        {
            score_index_vec.push_back(std::make_pair(scores[i], i));
        }
    }

    // Sort the score pair according to the scores in descending order
    std::stable_sort(score_index_vec.begin(), score_index_vec.end(),
        SortScorePairDescend<int>);

    // Keep top_k scores if needed.
    if (top_k > 0 && top_k < (int)score_index_vec.size())
    {
        score_index_vec.resize(top_k);
    }
}

template <typename BoxType>
inline void NMSFast_(const std::vector<BoxType>& bboxes,
    const std::vector<float>& scores, const float score_threshold,
    const float nms_threshold, const float eta, const int top_k,
    std::vector<int>& indices, float (*computeOverlap)(const BoxType&, const BoxType&))
{
    // Get top_k scores (with corresponding indices).
    std::vector<std::pair<float, int> > score_index_vec;
    GetMaxScoreIndex(scores, score_threshold, top_k, score_index_vec);

    // Do nms.
    float adaptive_threshold = nms_threshold;
    indices.clear();
    for (size_t i = 0; i < score_index_vec.size(); ++i) {
        const int idx = score_index_vec[i].second;
        bool keep = true;
        for (int k = 0; k < (int)indices.size() && keep; ++k) {
            const int kept_idx = indices[k];
            float overlap = computeOverlap(bboxes[idx], bboxes[kept_idx]);
            keep = overlap <= adaptive_threshold;
        }
        if (keep)
            indices.push_back(idx);
        if (keep && eta < 1 && adaptive_threshold > 0.5) {
            adaptive_threshold *= eta;
        }
    }
}


template <typename T>
static inline float rectOverlap(const T& a, const T& b)
{
    return 1.f - static_cast<float>(jaccardDistance<double>(a, b));
}

void NMSBoxes(const std::vector<Rectf>& bboxes, const std::vector<float>& scores,
    const float score_threshold, const float nms_threshold,
    std::vector<int>& indices, const float eta, const int top_k)
{
    NMSFast_(bboxes, scores, score_threshold, nms_threshold, eta, top_k, indices, rectOverlap);
}