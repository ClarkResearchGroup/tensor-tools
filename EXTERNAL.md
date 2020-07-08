Some portions of this code are modified from the ITensor library, which has the original copyright

```C++
//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
```

The list of modified files that include such work are:
- `linalg/tensor_davidson.cpp`
  - `tensor_davidsonIT`
- `models`
  - `hams/AutoMPO.cpp`
  - `hams/Heisenberg.cpp`
    - `buildHam(AutoMPO& ampo, MPO<T>&  H)`
  - `lattice`
    - `latticebond.h`,`triangular.h`,`square.h`
