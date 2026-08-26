# Autoresearch adaptation for the MVA hackathon

This is a constrained adaptation of the Karpathy autoresearch pattern. The
agent may propose changes to a ranking configuration or an annotation parser;
the evaluation harness is fixed and the challenge's hidden clinical answer is
never placed in the repo or used as a target.

Each experiment must:

1. run the unit tests and the fixed proxy checks;
2. record the exact config, code revision, elapsed time, and metrics;
3. explain why the change should generalise to a single rare-disease case;
4. keep a change only if it improves the proxy aggregate without regressing
   deterministic safety checks; and
5. stop after a small, predeclared experiment budget.

Live leaderboard submissions are manual, capped at six, and are not part of an
overnight autonomous loop. A leaderboard score is an external observation, not
a license to tune to the confirmed answer or to include private data in Git.

