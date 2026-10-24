// Hierarchical clustering of conformations (no DOM; also used by tests).

// Average-linkage (UPGMA) clustering: leaf order and dendrogram merges [a, b, height, size].
export function upgma(D) {
  const n = D.length;
  let clusters = Array.from({ length: n }, (_, i) => ({ id: i, members: [i], h: 0 }));
  const merges = [];
  const dist = (a, b) => { let s = 0; for (const i of a.members) for (const j of b.members) s += D[i][j]; return s / (a.members.length * b.members.length); };
  let next = n;
  while (clusters.length > 1) {
    let best = null;
    for (let i = 0; i < clusters.length; i++) for (let j = i + 1; j < clusters.length; j++) {
      const d = dist(clusters[i], clusters[j]);
      if (!best || d < best.d) best = { i, j, d };
    }
    const a = clusters[best.i], b = clusters[best.j];
    const c = { id: next++, members: [...a.members, ...b.members], h: best.d, left: a, right: b };
    merges.push(c);
    clusters = clusters.filter((_, k) => k !== best.i && k !== best.j);
    clusters.push(c);
  }
  return { root: clusters[0], order: clusters[0].members };
}

