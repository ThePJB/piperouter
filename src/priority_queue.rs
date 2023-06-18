use std::collections::HashMap;
use std::hash::Hash;

// todo remove index map maybe faster

pub struct PriorityQueue<K: Ord + Clone, V: Clone + Eq + Hash,> {
    pub heap: Vec<(K, V)>,
    pub index_map: HashMap<V, usize>,
}

impl <K: Ord + Clone,V: Clone + Eq + Hash> PriorityQueue<K, V> {
    pub fn new() -> PriorityQueue<K, V> {
        PriorityQueue { heap: Vec::new(), index_map: HashMap::new() }
    }

    pub fn push(&mut self, k: K, v: V) {
        if self.index_map.get(&v).is_none() {
            self.index_map.insert(v.clone(), self.heap.len());
            self.heap.push((k, v));
            self.upheap(self.heap.len() - 1);
            return;
        }
        let idx = *self.index_map.get(&v).unwrap();
        let (old_k, _) = self.heap[idx].clone();
        if k < old_k {
            self.heap[idx] = (k, v);
            self.upheap(idx);
        }
    }

    pub fn pop(&mut self) -> Option<(K, V)> {
        if self.heap.len() == 0 {
            return None;
        }
        let return_val = self.heap[0].clone();
        self.index_map.remove(&return_val.1);
        self.heap[0] = self.heap[self.heap.len() - 1].clone();
        self.index_map.insert(self.heap[0].1.clone(), 0);
        self.heap.truncate(self.heap.len() - 1);
        self.downheap(0);
        Some(return_val)
    }

    fn downheap(&mut self, mut idx: usize) {
        loop {
            let l = idx * 2 + 1;
            let r = idx * 2 + 2;

            if r < self.heap.len() {
                if self.heap[l].0 < self.heap[r].0 {
                    if self.heap[l].0 < self.heap[idx].0 {
                        self.swap(l, idx);
                        idx = l;
                    } else {
                        break;
                    }
                } else {
                    if self.heap[r].0 < self.heap[idx].0 {
                        self.swap(r, idx);
                        idx = r;
                    } else {
                        break;
                    }
                }
            } else if l < self.heap.len() {
                if self.heap[l].0 < self.heap[idx].0 {
                    self.swap(l, idx);
                    idx = l;
                } else {
                    break;
                }
            } else {
                break;
            }
        }
    }
    
    fn upheap(&mut self, mut idx: usize) {
        let parent = idx / 2;
        while self.heap[parent].0 < self.heap[idx].0 {
            self.swap(idx, parent);
            if parent == 0 {
                break;
            }
            idx = parent;
        }
    }

    fn swap(&mut self, i: usize, j: usize) {
        self.heap.swap(i, j);
        self.index_map.insert(self.heap[j].1.clone(), j);
        self.index_map.insert(self.heap[i].1.clone(), i);
    }
}