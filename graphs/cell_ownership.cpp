#include "cell_ownership.hpp"
#include <algorithm>
#include <unordered_set>

CellIndex AncestorCell(const CellIndex &cell, int current_level, int coarse_level)
{
    // Each coarsening halves the cell index, so the ancestor is
    // floor(i / 2^(current_level - coarse_level)) componentwise. C++ integer
    // division truncates toward zero, so negative indices are rounded down
    // explicitly.
    const int div = 1 << (current_level - coarse_level);
    CellIndex ancestor;
    for (int d = 0; d < 3; ++d)
    {
        const int i = cell[d];
        ancestor[d] = (i >= 0) ? i / div : -((-i + div - 1) / div);
    }
    return ancestor;
}

int ChildIndex(const CellIndex &cell)
{
    // Lexicographic child position (x fastest) to refinement ordering. The
    // first 4 entries are the quad ordering; all 8 are the hex ordering.
    static constexpr int map_from_lex[8] = {0, 1, 3, 2, 4, 5, 7, 6};
    const CellIndex parent = AncestorCell(cell, 1, 0);
    const int dx = cell[0] - 2 * parent[0];
    const int dy = cell[1] - 2 * parent[1];
    const int dz = cell[2] - 2 * parent[2];
    return map_from_lex[dx + 2 * dy + 4 * dz];
}

CellOwnership::CellOwnership(MPI_Comm comm_, const vector<CellIndex> &fine_cells,
                             int nlevels_)
    : comm(comm_), nlevels(nlevels_), owner(max(nlevels_ - 1, 0)),
      contributors(max(nlevels_ - 1, 0))
{
    int nranks, rank;
    MPI_Comm_size(comm, &nranks);
    MPI_Comm_rank(comm, &rank);

    const int fine_level = nlevels - 1;
    // One coarse level per coarsening: 0, ..., fine_level - 1
    const int ncoarsenings = fine_level;
    if (ncoarsenings <= 0) { return; }

    // Directory rank of a coarse cell. For now the root collects every cell;
    // a distributed directory only needs to change this function.
    auto directory_rank = [](int level, const CellIndex &cell) { return 0; };

    // 1. Send each distinct (level, ancestor cell) to its directory rank.
    // Sends vectors [level, i, j, k]
    vector<vector<int>> requests(nranks);
    {
        vector<unordered_set<CellIndex>> seen(ncoarsenings);
        for (const CellIndex &fine_cell : fine_cells)
        {
            // First coarse level to the coarsest.
            for (int level = ncoarsenings - 1; level >= 0; --level)
            {
                const CellIndex cell = AncestorCell(fine_cell, fine_level, level);
                // If this ancestor was already seen, so were all coarser ones.
                if (!seen[level].insert(cell).second) { break; }
                vector<int> &buf = requests[directory_rank(level, cell)];
                buf.insert(buf.end(), {level, cell[0], cell[1], cell[2]});
            }
        }
    }
    vector<int> request_counts;
    const vector<int> recv_requests = Exchange(comm, requests, request_counts);
    requests.clear();

    // 2. Directory assigns owner to each coarse cell. All contributing ranks
    //    are collected in `directory[level][cell]`. The first is the
    //    owner because contributing ranks are in ascending order.
    vector<unordered_map<CellIndex, vector<int>>> directory(ncoarsenings);
    {
        int pos = 0;
        for (int src = 0; src < nranks; ++src)
        {
            const int end = pos + request_counts[src];
            for (; pos < end; pos += 4)
            {
                const int level = recv_requests[pos];
                const CellIndex cell(recv_requests[pos + 1], recv_requests[pos + 2],
                                  recv_requests[pos + 3]);
                directory[level][cell].push_back(src);
            }
        }
    }

    // 3. Reply to each contributor with the owner; the owner also receives the
    // contributor list.
    // Sends [level, i, j, k, owner, n, contributor_1, ..., contributor_n]
    // with n = 0 for non-owners.
    vector<vector<int>> replies(nranks);
    for (int level = 0; level < ncoarsenings; ++level)
    {
        for (const auto &[cell, ranks] : directory[level])
        {
            const int cell_owner = ranks.front();
            for (int r : ranks)
            {
                vector<int> &buf = replies[r];
                const int n = (r == cell_owner) ? (int)ranks.size() : 0;
                buf.insert(buf.end(), {level, cell[0], cell[1], cell[2], cell_owner, n});
                if (n > 0) { buf.insert(buf.end(), ranks.begin(), ranks.end()); }
            }
        }
    }
    directory.clear();

    vector<int> reply_counts;
    const vector<int> recv_replies = Exchange(comm, replies, reply_counts);
    replies.clear();

    for (size_t pos = 0; pos < recv_replies.size(); )
    {
        const int level = recv_replies[pos];
        const CellIndex cell(recv_replies[pos + 1], recv_replies[pos + 2],
                          recv_replies[pos + 3]);
        const int cell_owner = recv_replies[pos + 4];
        const int n = recv_replies[pos + 5];
        pos += 6;

        owner[level][cell] = cell_owner;
        if (n > 0)
        {
            MFEM_VERIFY(cell_owner == rank, "contributor list sent to non-owner");
            contributors[level][cell].assign(recv_replies.begin() + pos,
                                          recv_replies.begin() + pos + n);
            pos += n;
        }
    }
}
