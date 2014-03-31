#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Author: Guillaume Filion.
Largely inspired from a blog post by Nick Johnson.
http://blog.notdot.net/2010/07/Damn-Cool-Algorithms-Levenshtein-Automata
Original gist code: https://gist.github.com/Arachnid/491973
"""

class NFA(object):
   """Non deterministic Finite state Automaton.

   The automaton is non deterministic in the sense that the same
   input symbol can cause several different states transitions.
   The special input symbol epsilon corresponds to transitions
   that do not consume any input symbol.

   The transitions map a source state and an input symbol to a
   set of possible destination states. Transitions are implemented
   as a nested dictionary. The first level associates a source state
   to a dictionary representing all possible transitions. The second
   level associates an input symbol to a set of destination states.

   The states are type agnostic as long as they are different
   Python objects."""

   # Special input symbols.
   EPSILON = object()
   ANY = object()
  
   def __init__(self, start_state):
      """Create an NFA with start state and no transition."""
      self.transitions = dict()
      self.final_states = set()
      self.start_state = start_state

   @property
   def states(self):
      """Output the set of all states of the NFA."""
      return set(self.transitions.keys())
 
   def add_transition(self, src, c, dest):
      """Add given transition to the NFA. Since the automaton
      is non deterministic, destinations are always implemented
      as a set of states."""
      self.transitions.setdefault(src, {}).setdefault(c, set()).add(dest)

   def add_final_state(self, state):
      """Specify that given state is final (accepting)."""
      self.final_states.add(state)
 
   def get_dest_states(self, src, c):
      """Return all possible transitions from state and symbol."""
      return self.transitions.get(src, {}).get(c, set())
  
   def to_dfa(self):
      """Return the powerset construction of the NFA."""

      # Helper functions.
      def plus_epsilon(stateset):
         """Extend state set by adding all states reachable by
         epsilon transitions."""
         extended_stateset = set(stateset)
         frontier = set(stateset)
         while frontier:
            state = frontier.pop()
            epsilon_moves = self.get_dest_states(state, NFA.EPSILON)
            # Add state only if has not been visited.
            new_stateset = epsilon_moves.difference(extended_stateset)
            frontier.update(new_stateset)
            extended_stateset.update(new_stateset)
         return extended_stateset

      def get_dest_stateset(src_stateset, c):
         """Return the destination state set from source state set
         and input symbol."""
         dest_stateset = set()
         for state in src_stateset:
            dest_stateset.update(self.get_dest_states(state, c))
            dest_stateset.update(self.get_dest_states(state, NFA.ANY))
         return frozenset(plus_epsilon(dest_stateset))

      def valid_symbols_from(stateset):
         """Return all valid symbols (to the possible excluding of
         epsilon) from given state set."""
         inputs = set()
         for state in stateset:
            inputs.update(self.transitions.get(state, {}).keys())
         # Exclude epsilon from valid input symbols.
         return inputs.difference([NFA.EPSILON])
 
      def add_transition_to(dfa, source, c, dest):
         """Add transition or default transition to DFA."""
         if c is NFA.ANY:
            dfa.set_default_transition(source, dest)
         else:
            dfa.add_transition(source, c, dest)

      def is_final(stateset):
         """Return True if any of the states in set is final."""
         intersection = self.final_states.intersection(stateset) 
         return True if intersection else False


      # Powerset construction. Starting from the initial state
      # set, determine all reachable state sets by merging all
      # states reachable from a state set and a given input
      # symbol. At each step, the state set is extended by the
      # epsilon transitions.
      start_stateset = frozenset(plus_epsilon(set([self.start_state])))
      dfa = DFA(start_stateset)

      seen = set()
      frontier = set([start_stateset])
      while frontier:
         stateset = frontier.pop()
         # Merge all states reachable from state set through
         # given input symbol.
         for c in valid_symbols_from(stateset):
            new_stateset = get_dest_stateset(stateset, c)
            add_transition_to(dfa, stateset, c, new_stateset)
            if new_stateset in seen: continue
            # Update the record if state set is seen for the
            # first time.
            frontier.add(new_stateset)
            seen.add(new_stateset)
            if is_final(new_stateset):
               dfa.add_final_state(new_stateset)

      return dfa


class DFA(object):
   """Deterministic Finite state Automaton.

   The automaton is deterministic in the sense that a given input
   symbol always causes only one state transitions. Every transition
   consumes an input symbol (there is no epsilon symbol).

   As for NFA, transitions are implemented as a nested dictionary.
   The difference is that the second level associates an input symbol
   to a single destination state. DFA also contain default transitions
   between states. The DFA will do a default state transition when
   the input symbol is valid but absent from the transitions. 

   As for NFA, the states are type agnostic."""

   def __init__(self, start_state):
      """Create an DFA with start state, no transition and no
      default transitions."""
      self.start_state = start_state
      self.transitions = {}
      self.defaults = {}
      self.final_states = set()
  
   def add_transition(self, src, c, dest):
      """Add transition to the DFA, overwritng when applicable."""
      self.transitions.setdefault(src, {})[c] = dest
  
   def set_default_transition(self, src, dest):
      """Add default transition, , overwritng when applicable."""
      self.defaults[src] = dest
  
   def add_final_state(self, state):
      """Specify that given state is final (accepting)."""
      self.final_states.add(state)

   def is_final(self, state):
      """Return whether state is final (accepting)."""
      return state in self.final_states
  
   def get_dest_state(self, src, c):
      """Return the destination state set from source state set
      and input symbol."""
      state_transitions = self.transitions.get(src, {})
      return state_transitions.get(c, self.defaults.get(src, None))


class LevenshteinAutomaton(object):
   """Levenshtein automaton. The states are represented as a pair
   of integers representing the coordinates in the match matrix."""

   def __init__(self, term, maxtau, substring=False):
      """Create an Levenshtein automaton from a term and a maximum
      edit distance maxtau. Also specify whether the automaton
      should match full or partial strings."""
      self.term = term
      self.matxau = maxtau
      self.substring = substring

      # Create the NFA and then the DFA from term and maxtau.
      # This will be the 'automaton' attribute of the object.
      nfa = NFA((0,0))
      # If matching in substring mode, add a universal transition
      # from the initial state to itself.
      if substring: nfa.add_transition((0,0), NFA.ANY, (0,0))

      for (i,c) in enumerate(term):
         nfa.add_transition((i,maxtau), c, (i+1,maxtau))
         for tau in range(maxtau):
            # Correct character
            nfa.add_transition((i,tau), c, (i+1,tau))
            # Substitution
            nfa.add_transition((i,tau), NFA.ANY, (i+1, tau+1))
            # Deletion
            nfa.add_transition((i,tau), NFA.ANY, (i, tau+1))
            # Insertion
            nfa.add_transition((i,tau), NFA.EPSILON, (i+1, tau+1))
               
      # The length of the term is equal to i+1.
      nfa.add_final_state((i+1,maxtau))
      for tau in range(maxtau):
         nfa.add_transition((i+1,tau), NFA.ANY, (i+1,tau+1))
         nfa.add_final_state((i+1,tau))

      self.automaton = nfa.to_dfa()


   def find(self, string):
      """Return the end position of the first match in the string,
      or -1 if no match is found."""
      align = []
      state = self.automaton.start_state
      for i, x in enumerate(string):
         state = self.automaton.get_dest_state(state, x)
         if self.automaton.is_final(state): return i
         if not state: break
      return -1

class PatternMatcher(object):
   def __init__(self, pattern, maxtau):
      self.pattern = pattern
      self.maxtau = maxtau
      self.forward_partial = LevenshteinAutomaton(pattern, maxtau, True)
      self.reverse_partial = LevenshteinAutomaton(reversed(pattern), maxtau, True)
      self.forward = LevenshteinAutomaton(pattern, maxtau)
      self.reverse = LevenshteinAutomaton(reversed(pattern), maxtau)

   def end(self, seq):
      start = self.reverse_partial.find(seq[::-1])
      return len(seq)-start-1 + self.forward.find(seq[len(seq)-start-1:]) if start > -1 else -1

   def start(self, seq):
      end = self.forward_partial.find(seq)
      return end - self.reverse.find(seq[end::-1]) if end > -1 else -1
