# 다음 실행의 재현성 검토 — 2026-09-17

판정: **StructCooker 전체가 다음 실행에서도 안정적으로 재현된다고 승인하기 어렵다.**
이번 데이터 릴리스의 진행률과 별도로, 정식 실행 경로·중단 후 재개·입력 정체성·결과
판정에 문제가 있다. PDB 전용 `pdb-build`에는 더 강한 보호가 있으나 일반 `build-all`과
현재 distillation 운영 스크립트에는 같은 보장이 적용되지 않는다.

현재 소스와 운영 스크립트를 읽고, 격리된 SLURM 잡 **268394**에서 작은 임시 입력과
mock으로 다섯 가지 동작을 재현했다. 실제 분자 데이터나 실행 중인 빌드 잡은 변경하지
않았다. 증거: `logs/reproducibility_review/20260917/{probe.py,findings.json,probe_268394.out}`.
전체 신규 데이터셋 재빌드나 새 환경 설치 실험은 수행하지 않았다.

## P1 — 완료 여부를 index 존재로 판단

근거: `libs/datacooker/src/datacooker/conditions.py:34`,
`src/structcooker/cli.py:443`.

일반 `build-all`은 대상의 `.index.tsv` 존재를 완료 신호로 사용한다. DB 자체가 없는
상태에서 stale index 파일만 생성해도 `output_absent=False`가 반환됐다. 입력·레시피가
바뀌었거나 미완성 산출물이 남은 경우에도 이 경로에는 검증된 receipt 확인이 없다.
Projection은 출력 파일 존재만 확인한다.

필요한 완료 조건: 입력/코드/설정 identity, 산출물 fingerprint, 필수 검증 receipt가 모두
맞을 때만 재사용해야 한다. `pdb-build`의 `BuildState` 보호를 일반 경로와 혼동하면 안 된다.

## P1 — distillation 전체 실행과 재개가 일회성 운영 스크립트에 의존

근거: `logs/distillation_followup/20260915/submit.py:6`, `work.py:11`,
`work.py:117`, `publish.py:29`, `setup.py:18`, `.gitignore:1`.

실행에 이전 잡 ID, 20260911/14 산출물, 고정 경로, 14개 기존 long shard와 18개 새 partition이
필요하다. Short/RNA 기본 MSA도 기존 파일을 전제로 한다. 따라서 새 출력 디렉터리에서
전체를 생성하는 동일한 실행 경로가 아니다. 필수 스크립트와 복구 참조가 `logs/` 아래에
있어 일반적인 checkout으로 전달되지 않는다. 라이브러리 소스 snapshot은 만들지만 실제
실행되는 `work.py`와 게시 스크립트는 그 snapshot에 포함되지 않는다.

게시 도중 DB 생성 후 index/receipt 작성 전에 중단되면, 재실행 시 `target.exists()`에서
거부되어 미완료 단계만 자동 재개하지 못한다. `SC_RESUME_BUILT`도 기존 build-report의
건수 일치만 먼저 확인하며, 당시 원본·설정·레시피 identity를 대조하지 않는다. 뒤의 payload
검증만으로 같은 key/shape를 유지한 다른 입력 버전까지 식별할 수는 없다.

필요한 완료 조건: 입력 release와 run directory를 받는 정식 native SLURM DAG, 새 빌드와
재개의 동일 경로, runner/config/reference까지 포함한 실행 snapshot, 단계별 receipt 및
게시 중단 지점별 복구. 상주 컨트롤러가 필요한 문제는 아니다.

## P1 — 숨은 메타데이터 입력과 프로세스 캐시

근거: `src/structcooker/instructions/readers/openfold.py:42`, `:57`, `:90`,
`src/structcooker/preflight.py:74`, `src/structcooker/production.py:168`.

서열-ID map과 구조 복구 manifest가 환경변수를 통해 reader 내부에서 선택되지만 해당
경로는 현재 YAML 입력 선언에 포함되지 않는다. 실제 `long_msa.yaml`에 대한
`preflight.input_paths` 결과는 raw data directory뿐이었다. 이 경로에만 의존하는 사전
검사와 fingerprint에는 별도의 ID map/복구 파일 정체성이 들어가지 않는다.

같은 프로세스에서 `MONOMER_SEQID_MAP`을 P1 map에서 P2 map으로 바꾼 뒤 같은 entry의 key를
다시 계산했는데 P1이 반환됐다. 전역 캐시가 경로나 내용 버전으로 구분되지 않는다.
현재 worker마다 새 프로세스를 쓰는 경우에는 이 캐시 증상이 나타나지 않을 수 있지만,
동일 프로세스의 다중 실행/API 사용에는 잘못된 key를 줄 수 있다.

필요한 완료 조건: ID map, 복구 manifest와 복구 payload를 명시적 설정 입력으로 선언하고
fingerprint에 포함하며, 캐시를 입력 identity로 구분해야 한다.

## P1 — 중복 MSA key 처리 결과가 완료 순서에 의존

근거: `libs/datacooker/src/datacooker/lmdb/core.py:993`,
`libs/datacooker/src/datacooker/_ray.py:163`.

공통 write 경로는 `txn.put(key, payload)`로 기존 값을 덮어쓴다. 동일 key의 A/B payload를
순서만 바꿔 전달한 검증에서 저장값이 B/A로 달라졌다. 양쪽 모두 written=2, DB entries=1,
failed=0이었다. Ray 완료 순서는 보장된 입력 순서가 아니므로 같은 서열의 입력이 둘 이상인
정식 build 경로는 스케줄링에 따라 선택되는 MSA가 달라질 수 있다.

현재 날짜별 distillation inventory는 별도 dedup으로 이 문제를 피한다. 따라서 이번
산출물이 전부 비결정적이라는 주장은 아니다. 다만 그 정책이 공통 build 계약으로 올라오지
않아 다음 실행 경로에 따라 보장이 달라진다.

필요한 완료 조건: 중복을 명시적으로 거부하거나 안정적인 canonical source 선택/병합
정책을 실행 전에 적용하고, ledger에서 중복과 최종 고유 key 수를 대조해야 한다.

## P1 — template 변환 오류가 정상적인 빈 결과로 통과

근거: `src/structcooker/instructions/transforms/template.py:969`.

chain key 부재는 `continue`, 디코딩/변환 예외도 광범위한 `except Exception: continue`로
처리된다. 손상된 chain payload를 모사했을 때 예외가 전파되거나 실패 목록이 생성되지
않고 `{}`를 반환했다. 이 결과를 감싼 `template_mols` 레코드는 schema D 검증을 통과했다.
정상적으로 hit가 없는 경우와 변환 실패로 모든 hit를 잃은 경우를 구분할 수 없다.

필요한 완료 조건: query별 후보 수·성공 수·제외 사유·원본 실패를 ledger로 남기고,
예상 가능한 제외와 처리 오류를 구분해야 한다. 완료 검증이 record 수와 shape만으로
이 차이를 놓치지 않아야 한다.

## P2 — 공유 노드 메모리 상태에 따라 동시성이 1에 머무를 수 있음

근거: `libs/datacooker/src/datacooker/_ray.py:36`, `:49`, `:133`, `:154`.

노드 전체 메모리 사용률이 60% 미만이어야 동시성이 증가한다. 62%에서 작업 완료를
10,000회 반영한 격리 검증에서도 최대 32 CPU 설정의 target은 1이었다. 65% 이상이면
추가 작업 admission도 막고 하나만 진행한다. Ray object store 크기를 명시하지 않고
각 잡이 개별 Ray를 시작하는 점도 공유 노드에서 자원 간섭을 키울 수 있다.

이는 동일한 입력에서 완료 시간과 timeout 성공 여부가 다른 잡 배치에 크게 좌우될 수
있음을 뜻한다. 실제 병목의 모든 원인을 이 단일 probe로 입증한 것은 아니다.

필요한 완료 조건: 잡별 메모리 예산, object-store/scratch 상한, 공유 노드 배치 정책,
가용 메모리에 맞춘 동시성 회복과 과부하 테스트. 단순히 요청 CPU만 늘리는 것으로
해결되지 않는다.

## P2 — 새 환경의 의존성 버전 재현이 repository만으로 고정되지 않음

근거: `.gitignore:6`, `pyproject.toml:9`, `src/structcooker/production.py:66`.

`pixi.lock`은 ignore 대상이고 현재 tracked 파일이 아니다. 여러 핵심 의존성이 범위 또는
와일드카드로 선언되어 있다. 실행 snapshot은 현재 lock과 설치 버전 정보를 기록하지만,
그 파일을 새 환경으로 전달·복원하는 공식 release 경로까지 보장하지는 않는다.

필요한 완료 조건: 배포에 lock 또는 동등한 고정 환경 artifact를 포함하고, 새 환경에서
설치와 작은 전체 DAG 재현을 확인해야 한다. 현재 작업 트리의 미커밋 상태 자체를
소프트웨어 결함으로 센 것은 아니다.

## 이미 있는 보호와 검토 한계

- `pdb-build`에는 입력/산출물 content fingerprint, 상태 receipt, 출력 소유권과 잠금,
  검증된 결과 재사용, 불명확한 제출 상태에서 중복 제출을 막는 처리가 있다.
- 5M 기준 자체와 shard 간 중복 key 검사는 구현되어 있다. 이번 문제를 5M 정책 탓으로
  볼 근거는 없다.
- 현재 177개 테스트 통과는 위의 운영 조건을 모두 포괄했다는 의미가 아니다. 이번
  검토는 소규모 동작 재현과 정적 분석이며 장시간 timeout, 노드 장애, 동시 게시 중단을
  실제 production 데이터에 주입하지 않았다.
- 재현 목표는 같은 입력/설정에 대한 key·feature·제외 결과의 일치여야 한다. 현재 codec은
  array 식별자로 UUID를 사용하므로 독립 실행의 압축 byte SHA가 같다는 조건과 의미상
  결과 일치를 혼동하면 안 된다.

승인 순서: 잘못된 성공 판정과 데이터 손실부터 차단하고, 입력 identity와 공식 재개
경로를 통합한 뒤, 작은 고정 fixture의 새 빌드/재실행/중단 복구 및 공유 노드 실행을
검증한다. 이 조건이 충족되기 전에는 이번 릴리스의 성공만으로 도구 전체를 프로덕션
수준이라고 선언하지 않는다.

## 2026-09-17 수정 및 검증

앞의 항목은 수정 전 감사 결과이다. 후속 수정은 다음과 같다.

- `build-all`을 `release.py`의 유한 SLURM 단계 제출 경로로 전환했다. DB/index
  존재만으로 완료를 판단하지 않고 입력·코드·출력 identity와 전체 LMDB 검증 receipt를
  사용한다. 동작과 재개 방법은 [release-build.md](release-build.md)에 기록했다.
- 입력 중복은 기본 거부한다. distillation MSA의 명시적 `first_path` 정책은 정렬된
  경로로 결정하며 제외 경로를 기록한다. 엔진의 중복 output key도 오류로 처리한다.
- reader map/override를 config 입력으로 선언하고 attempt에 복사한다. 직접 호출의
  process cache도 파일 경로 및 stat identity 변화에 따라 갱신한다.
- template decode/transform 오류를 전파하고, 누락 chain과 선택 제외를 레코드 ledger에
  남긴다. schema D가 ledger의 후보 수와 key 집합 일치를 검사한다.
- admission 회복 임계값을 조정하고 SLURM 메모리 예산에 따른 Ray object store 크기를
  명시했다. 62% 부하에서 동시성 회복, 80%에서 감소하는 회귀 검증을 추가했다.
- `pixi.lock`을 ignore에서 제외하고 SLURM 실행에 `--frozen`을 적용했다. 릴리스 때
  루트 변경과 DataCooker 서브모듈 변경, lock을 함께 배포해야 한다.

실제 SLURM 소규모 검증: `logs/reproducibility_fix/20260917/acceptance-run-v2`.
첫 실행 268402–268411에서 원본 MSA → cap → 검증 완료, 268413–268414에서
attempt를 그대로 재사용했다. fixture cap index를 삭제한 뒤 268415–268419에서
기존 shard로 재게시·인덱스·검증을 완료했다. production DB에는 장애를 주입하지 않았다.
회귀 테스트 job 268412: **190 passed**. Ruff 전체 검사 통과, Pyright 오류/경고 0.

남은 검증 범위: 별도 새 환경의 lock 설치, 장시간/대규모 자원 경합 및 노드 장애
시험은 수행하지 않았다. projection 파일은 현재 비어 있지 않은지와 content identity를
검사하며 모든 파일 포맷의 의미 검증을 제공하지 않는다. 제출 결과가 불명확한 경우는
안전하게 중단하여 수동 확인을 요구한다. 이 한계까지 검증했다고 표현해서는 안 된다.

최종 MSA cap/원본 key 보존 검사 추가 후 job **268421**도 **190 passed**.
공식 distillation manifest의 26개 단계 dependency dry-run 및 최종 Ruff 검사 통과.
